#!/usr/bin/env Rscript
# # Automated to run from terminal, see:
#   * README_scr2.md *
# REQUIRES : results from 1_....R
# performs PCA on cpm, hierarchical clust on cpm 
# OUTPUT:  plots, and saves  objects into rds/ folder
# johaGL 2022
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(pheatmap)
library(patchwork)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
gv <- yaml::read_yaml(args[1])

outname <- gv$outname
equi_id <- gv$equivalentid  # inspect gene id's to know ('.'? --> _version) !
shortname <- gv$shortname

SPECIESensmbldsetname <- gv$SPECIESensmbldsetname
samplestodrop <- gv$samplestodrop
typeRNAs <- gv$typeRNAs # a list of biotypes
print(typeRNAs)
# START
setwd(gv$mywdir)
source(paste0(gv$mywdir,"scr2/func.R")) 
o_dir <- paste0(gv$mywdir,"results_",outname,"/")
resdirs <- getoutputdirs(outer=o_dir)
rds_dir = resdirs[1]; tabl_dir = resdirs[2]; plo_dir = resdirs[3] 


################### open readyData ################
readydafile <- paste0(rds_dir, "readyData",outname,".rds")

myobj <- readRDS(readydafile)

################### if typeRNAs then filter all object ################

if (is.null(typeRNAs) == TRUE){
  print("You have not defined typeRNAs, nothing to do")
}else{
  myobj <- o_getrows(myobj, 'biotype', typeRNAs )
}

#################### cpm ###################
myobj@cpm <- apply(myobj@rawq, 2, function(col) return((col/sum(col))*(1e+06)))
corner(myobj@cpm)

# # ----------------------------- Multivariate -----------------------
# factoextra and factominer
res.pca <- FactoMineR::PCA(t(myobj@cpm), scale.unit = FALSE, graph = FALSE)
saveRDS(res.pca, paste0(rds_dir, "cpm-pcaData_",outname,".rds") )

dists <- as.matrix(dist(t(myobj@cpm)))
onedists <- 1 - dists
corner(onedists)
condition_annotation <- myobj@metadata[,"condition"]
names(condition_annotation) <- myobj@metadata$sample
if (length(myobj@metadata$sample)>15){
  flexhierpdf <- round((length(myobj@metadata$sample)/2.5),0)
}else{flexhierpdf=6}
pdf(paste0(plo_dir,"cpm-hierarchical_", outname, ".pdf"), height=flexhierpdf+1,
    width=round(flexhierpdf*1.5),0)
pheatmap(onedists,
         main = paste("hierarchical clustering,", outname) ,
         annotation_col=as.data.frame(condition_annotation))
dev.off()
dev.off()


 
