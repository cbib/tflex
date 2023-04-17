#!/usr/bin/env Rscript
# # Automated to run from terminal, see:
#   * README_scr2.md *
# Initial verification, creates output directories,
# checks biotypes, and performs multivariate and univariate analyses
# OUTPUT: multiple plots, and saves dds object into rds/ folder
# usage:
#   Rscript $DIRTFLEX/Rscripts/2_a_DE.R config_postprocess.yml $DIRTFLEX/ 

# johaGL 2022
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(biomaRt)
library(DESeq2)
library(pheatmap)
library(patchwork)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
gv <- yaml::read_yaml(args[1])

conditionLEVELS <- gv$conditionLEVELS
reqcontrast <- gv$requiredcontrast 
objecthasthiscontrast <- gv$objecthasthiscontrast
outname <- gv$outname
equi_id <- gv$equivalentid  # inspect gene id's to know ('.'? --> _version) !
shortname <- gv$shortname

SPECIESensmbldsetname <- gv$SPECIESensmbldsetname
samplestodrop <- gv$samplestodrop

source(paste0(args[2], "Rscripts/func.R"))


# START
setwd(gv$mywdir)

o_dir <- paste0(gv$mywdir,"results_",outname,"/")
resdirs <- getoutputdirs(outer=o_dir)
rds_dir = resdirs[1]; tabl_dir = resdirs[2]; plo_dir = resdirs[3] 
strcontrast <- paste(reqcontrast[2:3],collapse="_vs_")


# ---------------------------- set auto files names  --------------------

print(objecthasthiscontrast)
if (is.null(objecthasthiscontrast)){  
  dds_finame <-  paste0(rds_dir, "ddsObj_",strcontrast,"_",outname,".rds")
  readydafile <- paste0(rds_dir, "readyData",outname,".rds")
  plotsmultirlog <- TRUE
}else{
  print("using dds object already saved")
  dds_finame <- objecthasthiscontrast
  readydafile <- paste0(rds_dir, "readyData",outname,".rds")
  plotsmultirlog <- FALSE
}
print(dds_finame); print(readydafile)

# ----------------------- open 'readyData' .rds object  ---------------------
myobj <- readRDS(readydafile)
if(all(rownames(myobj@rawq)==myobj@row_data$symbol_unique)){
  print("as readyData matrix has symbols, re-impute ensembl_id to this")
  rownames(myobj@rawq) <- myobj@row_data$x # see func.R : getfullRowsSymbols
}else{ print("conflict rownames and row_data in object, check 1_....R")}

# ---------------- Univariate simple contrast (DESeq2) -------------------
print("Performing only Univariate simple contrast, see README_scr2.md")
print("(LEVELS and nb libraries determine what you get: resultsNames(dds)")
cat('\n  * * * 1_initiotodds.R runs univariate analysis, or opens existing rds',
    '(I do not overwrite) * * * \n')

if (!file.exists(dds_finame )){
  print(head(rownames(myobj@rawq)))
  cat("\nset factors with levels and converting into DESeq2 object\n" )
  myobj@metadata$condition <- factor(myobj@metadata$condition, levels=conditionLEVELS)
  dds <- DESeqDataSetFromMatrix(countData = myobj@rawq,
                                colData = myobj@metadata,
                                design= ~ condition)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- DESeq(dds)
  saveRDS(dds,file=dds_finame )
}else{
  dds <- readRDS(dds_finame)
}
print("here your contrasts in this dds object")
print(resultsNames(dds))

# --------------- Estimates, and Multivariate with rlog matrix -----------------------
if (plotsmultirlog){
  # ------------------- Estimates ----------------------
  # pdf(paste0(plo_dir,"deseq2EstsMA_",strcontrast,"_",outname,".pdf"), 
  #     width=11, height=5)
  tiff(paste0(plo_dir,"deseq2EstsMA_",strcontrast,"_",outname,".tiff"), res = 80,
       width=12, height=6, units="in")
  par(mfrow = c(1,2))
  plotDispEsts(dds, main=paste("Dispersion Estimates,",outname))
  plotMA(dds, ylim=c(-5,5),  alpha=0.05, main=paste("MA plot,",outname) )
  par(mfrow=c(1,1))
  dev.off()
  # ----------------------- Multivariate with rlog matrix ----------------- 
  rloo <- rlog(dds)
  dists <- as.matrix(dist(t(assay(rloo))))
  onedists <- 1 - dists
  condition_annotation <- myobj@metadata[,"condition"]
  names(condition_annotation) <- myobj@metadata$sample
  if (length(myobj@metadata$sample)>15){
    flexhierpdf <- round((length(myobj@metadata$sample)/2.5),0)
  }else{flexhierpdf=6}
  pdf(paste0(plo_dir,"rlog-hier_", outname, ".pdf"), height=flexhierpdf+1,
      width=round(flexhierpdf*1.5),0)
  pheatmap(onedists, 
           main = paste("hierarchical clustering,", outname) ,
           annotation_col=as.data.frame(condition_annotation))
  dev.off()
  
  # factoextra and factominer
  res.pca <- FactoMineR::PCA(t(assay(rloo)), scale.unit = FALSE, graph = FALSE)
  saveRDS(res.pca, paste0(rds_dir, "rlog-pcaData_",outname,".rds") )
  screeplo <- fviz_eig(res.pca, addlabels=TRUE) + 
    labs(title = paste( "Scree plot,", outname ))
  
  pcaind <- fviz_pca_ind(res.pca, col.ind = "cos2", 
                         gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                         repel=TRUE, graph=FALSE)
  
  pdf(paste0(plo_dir,"rlog-pca_",outname,"_ind.pdf"), height=5, width=10)
  print(screeplo | pcaind  )#using patchwork
  dev.off()
}

# END


