#!/usr/bin/env Rscript
# # Automated to run from terminal, similar to 2_a_DE.R

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(patchwork)
library(RColorBrewer)
library(yaml)

args = commandArgs(trailingOnly=TRUE)
gv <- yaml::read_yaml(args[1])

reqcontrast <- gv$requiredcontrast 
objecthasthiscontrast <- gv$objecthasthiscontrast
outname <- gv$outname
equi_id <- gv$equivalentid
shortname <- gv$shortname
labelvolca <- gv$labelsinvolcano

if (is.null(labelvolca)){
  labelvolca <- TRUE
}

### vars  in 2_ddstotabplots.R
trustedcolumn <- gv$trustedcolumn  #useful: symbols, set 'symbol_unique' !
absLFCcutoff <- gv$absLFCcutoff 
sigcutoff <- gv$sigcutoff 
aLFClab <- gv$aLFClab 
pjlab <- gv$pjlab 

source(paste0(args[2], "Rscripts/func.R"))

# START
setwd(gv$mywdir)

o_dir <- paste0(gv$mywdir,"results_",outname,"/")
resdirs <- getoutputdirs(outer=o_dir)
rds_dir = resdirs[1]; tabl_dir = resdirs[2]; plo_dir = resdirs[3] 
strcontrast = paste(reqcontrast[2:3],collapse="_vs_")

if (trustedcolumn != 'symbol_unique'){
  print("error! define trustedcolumn as 'symbol_unique' in .yaml")
}

print(objecthasthiscontrast)
if (is.null(objecthasthiscontrast)){  
  dds_finame <-  paste0(rds_dir, "ddsObj_",strcontrast,"_",outname,".rds")
}else{
  dds_finame <- objecthasthiscontrast

}
readydafile <- paste0(rds_dir, "readyData",outname,".rds")

print(dds_finame); print(readydafile)

# ------------------- Open saved dds object rds, doing advanced plots--------------------
myobj <- readRDS(readydafile)
dds <- readRDS(dds_finame)
print(resultsNames(dds))
res.here <- DESeq2::results(dds, reqcontrast)
head(res.here)

# organize entire results 
exclude_types <- c("rRNA","snRNA")
print(paste("saving Extendend results BUT EXCLUDing : ",paste(exclude_types)))

print(dim(res.here))
print(all(res.here$ensembl_gene_id_version == unique(res.here$ensembl_gene_id_version)))

resExtended <- symbomore4DE(res.here, myobj, trustedcolumn, equi_id, 
                     exclude_types, tabl_dir, outname)

savetsv(resExtended,
            paste0(tabl_dir,"DEextended_",strcontrast,"_",outname,".tsv"))
# top features for saving to csv : 
print("saving top results [overwrites existing tables!]") 
rtop <- resExtended %>% filter(padj <= sigcutoff & 
                                 abs(log2FoldChange) >= absLFCcutoff) 
savetsv(rtop, paste0(tabl_dir,"DEtop_",strcontrast,"_",outname,".tsv"))

uu <- rtop %>% filter(log2FoldChange >= 0)
dodo <- rtop %>% filter(log2FoldChange < 0)
savetsv(uu, paste0(tabl_dir,"DEtop-UP_",strcontrast,"_",outname,".tsv"))
savetsv(dodo, paste0(tabl_dir,"DEtop-DOWN_",strcontrast,"_",outname,".tsv"))
# ------------------ VIZZZZZZ --------------------------------------
col_here <- c("royalblue","palevioletred3","firebrick") # cold --> hot color

print("preparing new dataframe for plot")
rall4plot <- symbolsasrownames(resExtended, trustedcolumn)
rall4plot <- reformat4plot_withlabels(rall4plot,
                                      rtop,
                                      sigcutoff,
                                      absLFCcutoff,
                                      trustedcolumn,
                                      aLFClab, pjlab, col_here)

if (labelvolca){
  pdf(paste0(plo_dir,"volcano_",strcontrast,"_",outname,".pdf"), height=16, width=8)
  print(cowplot::plot_grid(
    dovolcano_picklabcolumn(rall4plot, rall4plot$labelfew, strcontrast,
                            sigcutoff, absLFCcutoff,
                            pjlab, aLFClab,
                            outname,
                            mycaption="labeling for extreme padj or absLFC") ,
    dovolcano_picklabcolumn(rall4plot, rall4plot$labelrich, strcontrast,
                            sigcutoff, absLFCcutoff,
                            pjlab, aLFClab,
                            outname,
                            mycaption=paste("cutoffs for labels:",
                                            pjlab, ",",aLFClab,"; and extreme")),
    ncol = 1 ) # end plotgrid
  )#end print
  dev.off()
  }else{
    pdf(paste0(plo_dir,"volcano_",strcontrast,"_",outname,".pdf"))
    print(dovolcano_nolabel(rall4plot,  strcontrast,
                              sigcutoff, absLFCcutoff,
                              pjlab, aLFClab,  outname, mycaption="")
    )#end print
    dev.off()
  }

dev.off()


