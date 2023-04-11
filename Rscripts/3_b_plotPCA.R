#!/usr/bin/env Rscript
# # Automated to run from terminal

# performs PCA on cpm, hierarchical clustering plots 
# OUTPUT:  plots, and saves  objects into rds/ folder
# johaGL 2022
library(tidyverse)
library(FactoMineR)
library(factoextra)
#library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
#gv <- yaml::read_yaml("ratvars.yaml")
gv <- yaml::read_yaml(args[1])

conditionLEVELS <- gv$conditionLEVELS
# reqcontrast <- gv$requiredcontrast 
outname <- gv$outname
equi_id <- gv$equivalentid  # inspect gene id's to know ('.'? --> _version) !
shortname <- gv$shortname
SPECIESensmbldsetname <- gv$SPECIESensmbldsetname
samplestodrop <- gv$samplestodrop
cutoffs <- gv$cutoffsvarspvals
maxncontrib <- gv$maxncontrib # too extreme vars (genes) to exclude from biplot
dobiplots <- gv$dobiplots_bypval_bycond

source(paste0(args[2], "Rscripts/func.R"))

# START

#autopicker colors:

cond_colorpick <- function(shortname){
  print("twocond_colorpick : outputs a 8 length vector, different palettes by species")
  # colpicker <- list("mouse"=c("cadetblue","sienna2"), 
  #                   "rat"=c("lightblue4","salmon3"),
  #                    "human"=c("palegreen4","orange"))
  colpicker <- list("mouse"=brewer.pal(9, "Blues")[1:9],
                    "rat"=brewer.pal(9, "Greys")[1:19],
                    "human"=brewer.pal(9, "PuBuGn")[1:9])  
  CHOSENCOLORS <- colpicker[[shortname]] 
  return(rev(CHOSENCOLORS))
}

CHOSENCOLORS <- cond_colorpick(shortname)[1:length(conditionLEVELS)] 
names(CHOSENCOLORS) <- conditionLEVELS

setwd(gv$mywdir)

o_dir <- paste0(gv$mywdir,"results_",outname,"/")
resdirs <- getoutputdirs(outer=o_dir)
rds_dir = resdirs[1]; tabl_dir = resdirs[2]; plo_dir = resdirs[3] 

################### open readyData ################
readydafile <- paste0(rds_dir, "readyData",outname,".rds")

myobj <- readRDS(readydafile)


# # ----------------------------- open existing  PCA data-----------------------
res.pca <- readRDS( paste0(rds_dir, "cpm-pcaData_",outname,".rds") )

#res.pca <- readRDS("~/cocultureProj/results_Neurons_rat/rds/cpm-pcaData_Neurons_rat.rds")

screeplo <- fviz_eig(res.pca, addlabels=TRUE) + 
  labs(title = paste("Screeplot", outname))

twodimElipse <- fviz_pca_ind(res.pca, 
                              geom= c("point","text"), 
                              axes=c(1,2),
                        col.ind = as.factor(myobj@metadata$condition),
                         repel=TRUE,
                         addEllipses = TRUE,
                         invisible="quali",
                         palette= CHOSENCOLORS, 
                         graph=FALSE,
                         title = paste("PCA",outname) ) +
  theme(legend.position = "bottom")

twodimNoElipse <- fviz_pca_ind(res.pca, geom= c("point","text"), 
                               axes = c(1,2),
                               habillage = as.factor(myobj@metadata$condition),
                               repel=TRUE,
                               addEllipses = FALSE,
                               invisible="quali",
                               palette= CHOSENCOLORS,
                               graph=FALSE,
                               title = paste("PCA", outname)) + 
            theme(legend.position = "bottom")

boxplo <- doboxplot3PC_bicond(res.pca, myobj@metadata, 
                    CHOSENCOLORS,
          paste("Three first PCs (Dims) by condition",outname) )



pdf(paste0(plo_dir,"cpm-PCA-elip-",outname,".pdf"), height=5, width=7)
twodimElipse 
dev.off()

pdf(paste0(plo_dir,"cpm-PCA-noelip-",outname,".pdf"), height=5, width=7)
twodimNoElipse
dev.off()
 
pdf(paste0(plo_dir,"cpm-PCA-scree-",outname,".pdf"), height=5, width=5)
screeplo
dev.off()

pdf(paste0(plo_dir,"cpm-PCA-boxplot",outname, ".pdf"), height=4.5, width=6)
boxplo + theme(legend.position = "bottom")
dev.off()

# ------  PCA, explore vars  -----------------
coefis <- dimdesc(res.pca, axes=1)
coefdf <- as.data.frame(coefis$Dim.1$quanti)
coefdf$symbol_unique <- rownames(coefdf)

pdf(paste0(plo_dir, "cpm-PCAvarsDim1-hist_",outname,".pdf"))
ggplot(coefdf,aes(x=p.value)) + geom_histogram( fill="darkcyan",alpha=0.7) + theme_bw() +
  labs(title="p values Dim1 vars")
dev.off()
tooextreme <- names(rev(sort(res.pca$var$contrib[,"Dim.1"]))[1:maxncontrib])
df4plot <- coefdf %>% filter(!symbol_unique %in% tooextreme)

# ---- test cutoffs, i.e. do biplot, CD , viz purposes: NO too extreme --------------------
if (dobiplots == TRUE){
  lgg <- list()
  for (cutoff in cutoffs){    # tooextreme, cutoff,  meta, COLORS, df4plot
    exfiltp <- df4plot %>% filter(p.value <= cutoff)
    nameplot <- as.character(cutoff)
    lgg[[nameplot]] <- fviz_pca_biplot(res.pca,
                                       geom.ind="point",
                                       fill.ind=myobj@metadata$condition, #check
                                       col.ind="white",
                                       pointshape=21,
                                       pointsize=3,
                                       labelsize=2,
                                       palette=CHOSENCOLORS, #check
                                       addEllipses=FALSE,
                                       invisible="quali",
                                       repel=TRUE,
                                       col.var="grey20",
                                       alpha.var=0.6,
                                       #gradient.cols = "Greys",
                                       select.var=list(name=exfiltp$symbol_unique)) +
      theme(legend.position="bottom") +
      labs(fill="Condition",color="Contrib", alpha="Contrib",
           title = paste("PCA variables, p.value <=",cutoff, ":",dim(exfiltp)[1],"vars")
      ) 
    
  }
  
  
  annobar <- ggplot() + annotate("text",  x=1,y=1,size = 4,
                                 label = paste(" excluding",paste(tooextreme, collapse=", "))) + 
    theme_void()
  pdf(paste0(plo_dir,"cpm-PCAvars-",outname, ".pdf"),height = 5, width=5*length(lgg))
  plot_grid(plot_grid(plotlist = lgg, ncol=length(lgg)),
            annobar, nrow=2, rel_heights = c(5,0.5))
  dev.off()
  
}





