#!/usr/bin/env Rscript
# # Automated to run from terminal, see:
#   * README_scr2.md *
# Initial verification, creates output directories,
# checks biotypes,
# OUTPUT:  saves dds object into rds/ folder
# note : obj@row_data may have symbol_unique but empty external_gene_name: 
#                      because when symbol retrieved by ncbi , not biomart. 
# johaGL 2022
library(tidyverse)
library(biomaRt)
library(patchwork)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
gv <- yaml::read_yaml(args[1])

metadatafn <- gv$metadatafile
cfn <- gv$countsfile
conditionLEVELS <- gv$conditionLEVELS
reqcontrast <- gv$requiredcontrast 
objecthasthiscontrast <- gv$objecthasthiscontrast
outname <- gv$outname
equi_id <- gv$equivalentid  # inspect gene id's to know ('.'? --> _version) !
shortname <- gv$shortname

SPECIESensmbldsetname <- gv$SPECIESensmbldsetname
nb_bt <- gv$nb_bt # min nb of genes req for a biotype, for biotype to be displayed
samplestodrop <- gv$samplestodrop
typeRNAs <- gv$typeRNAs
# START
setwd(gv$mywdir)
source(paste0(gv$mywdir,"scr2/func.R")) 
o_dir <- paste0(gv$mywdir,"results_",outname,"/")
resdirs <- getoutputdirs(outer=o_dir)
rds_dir = resdirs[1]; tabl_dir = resdirs[2]; plo_dir = resdirs[3] 
strcontrast <- paste(reqcontrast[2:3],collapse="_vs_")

print("make sure samples names are ok, not being possible to change on the run")
print("if not the case, abort this and generate new matrix and metadata")
print("then put in your .yaml file those new paths ! ")
metadata <- read.table(metadatafn, header=TRUE,sep=',')
rawc <- opentsv(cfn) 
print_warnings_rawcdfANDmetadata(rawc,metadata)
print("")
############################ BIOTYPES LOAD    #############################
print(" * USING BIOMART FOR BIOTYPES AND MORE * ")
savebiomartinfos(SPECIESensmbldsetname, gv$mywdir, shortname)
bm_sp <- readRDS(paste0(gv$mywdir, "biomart_db/bm_",shortname,".rds"))
print("")
############## matrix cross with biomart , then  rnaexp object   ##############

rawmat <- countsdftomatrix(rawc)
rowdatadf <- getfullRowsSymbols(rawmat,equi_id, bm_sp) # precursor for @row_data

# ----------- further query, omit dot and digits after-----------------------------
if (equi_id=="ensembl_gene_id_version"){
  rowdatadf <- furtherRowsSymbols(rowdatadf, equi_id, bm_sp)
  write.csv(data.frame("searchingfurther_ensembl_gene_id"=savedesperate),
            paste0(tabl_dir,"querybiomartsuppl.csv"), row.names=FALSE)
  } else if (equi_id=="ensembl_gene_id"){
  print(paste(equi_id, "used ensembl_gene_id, no further query for rowdatadf"))
}

########### NEW : check symbols that were not retrieved, use ncbi ###############
nosym_db <- bm_sp %>% filter(is.na(external_gene_name) | external_gene_name == '')
# note: ncbi use only ensembl_gene_id (NEVER: _version)
opn <- save_ncbi_hits(nosym_db$ensembl_gene_id, shortname) # saves or get if done

nosymbol.here <- rowdatadf %>% 
  filter(is.na(external_gene_name) | external_gene_name == '')

ncbi.here <- opn %>% filter(ENSEMBL %in% nosymbol.here$ensembl_gene_id)
print(dim(ncbi.here))

K <- rowdatadf;   
if (shortname == "human"){
  tmp <- left_join(ncbi.here, nosym_db[ , c(equi_id,
                                    "ensembl_gene_id") ],
              by=c("ENSEMBL"="ensembl_gene_id") ) 
  ncbi.here <- tmp
  for (k in 1:dim(ncbi.here)[1]){
    v <- ncbi.here[k,]
    idi <- v[[equi_id]]
    
    if (dim(K[K$x==idi,])[1] > 0 & (!is.null(idi))) {
      K[K$x==idi, "symbol_unique"] <- v$SYMBOL
      if ( is.na(K[K$x==idi, "description"]) | (K[K$x==idi,"description"] == "") ){
        K[K$x==idi,"description"] <- v$GENENAME
      }
      if ( is.na(K[K$x==idi, "gene_biotype"]) | (K[K$x==idi,"gene_biotype"] == "") ){
        K[K$x==idi,"gene_biotype"] <- v$GENETYPE
      }
    } # end if dim K & null idi
  }
} else if(shortname == "rat") {
  for (k in 1:dim(ncbi.here)[1]){
    v <- ncbi.here[k,]
    idi <- v$ENSEMBL # not equi_id
    
    if (dim(K[K$x==idi,])[1] > 0 & (!is.null(idi))) {
      K[K$x==idi, "symbol_unique"] <- v$SYMBOL
      if ( is.na(K[K$x==idi, "description"]) | (K[K$x==idi,"description"] == "") ){
        K[K$x==idi,"description"] <- v$GENENAME
      }
      if ( is.na(K[K$x==idi, "gene_biotype"]) | (K[K$x==idi,"gene_biotype"] == "") ){
        K[K$x==idi,"gene_biotype"] <- v$GENETYPE
      }
    } # end if dim K & null idi
  }
}
K$symbol_unique <- make.unique(K$symbol_unique)
rowdatadf <- K

print(paste("rawmat is matrix",is.matrix(rawmat), "  and is array",  is.array(rawmat) ))
if(all(rownames(rawmat) == rowdatadf$x)){
  print("Creating object S4 type ('rnaexp' class), rownames match")
} else {
  print("WARNING: rownames wont match for object S4 type ('rnaexp' class)!!")
}
symat <- rawmat; rownames(symat) <- rowdatadf$symbol_unique
myobj <- new("rnaexp",rawq=symat, metadata=metadata, row_data=rowdatadf)
rm(symat,rawmat, rawc)
print("")

print("")



# --------------- keep/drop samples if said ----------------
if(is.null(samplestodrop)==FALSE){
  myobj <- o_getcols(myobj, "-", samplestodrop)
}

# ----------------- keep genes if said -----------------------:
if (is.null(typeRNAs) == TRUE){
  print("You have not defined typeRNAs (protein_coding, or other), nothing to do")
}else{
  myobj <- o_getrows(myobj, 'biotype', typeRNAs )
}
print("")
# -------------------  ** save object ** --------------------
print("saving object (class 'rnaexp') for future purposes (symbols inside)")
saveRDS(myobj,
         paste0(rds_dir, "readyData",outname,".rds") )


##################### plot biomart and OUR data  ####################### 
colnames(myobj@row_data)
titleX <- paste("Biotypes in whole", outname, "dataset")
X <- do_biotype_plot(myobj@row_data$gene_biotype, nb_bt, titleX )

pdf(paste0(plo_dir,"Biotypes_full", outname, ".pdf"), height= 4, width=6)
X 
dev.off()






   

