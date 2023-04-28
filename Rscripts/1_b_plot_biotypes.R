# Do the plot for the RNA biological type (protein, lnc,etc)
# usage:
#  Rscript $TFLEX/Rscripts/1_b_plot_biotypes.R config_postprocess.yml $TFLEX/


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

source(paste0(args[2], "Rscripts/func.R"))

# START
setwd(gv$mywdir)

o_dir <- paste0(gv$mywdir,"results_",outname,"/")
resdirs <- getoutputdirs(outer=o_dir)
rds_dir = resdirs[1]; tabl_dir = resdirs[2]; plo_dir = resdirs[3] 

do_biotype_plot <- function(gene_biotype, nb_bt, mytitle){

  biotyfreq <- table(gene_biotype)
  sortedfq <- sort(biotyfreq[biotyfreq >= nb_bt])
  bty <- data.frame("Biotype"=factor(names(sortedfq), levels=names(sortedfq)),
                    "Count"=as.vector(sortedfq)) 
  
 
  biotypes_colors <- colorRampPalette(brewer.pal(8,"Dark2"))(nrow(bty))

  gp <- ggplot(bty, aes(x=Biotype, y=Count, fill=biotypes_colors)) +
    geom_bar(stat="identity") +
    coord_flip() + theme_light() + theme(legend.position = "none") +
       labs(title = mytitle)
    
  return(gp)
}


##################### plot biomart and OUR data  ####################### 
readydafile <- paste0(rds_dir, "readyData",outname,".rds")
myobj <- readRDS(readydafile)

colnames(myobj@row_data)
titleX <- paste("Biotypes in whole", outname, "dataset")
X <- do_biotype_plot(myobj@row_data$gene_biotype, nb_bt, titleX )

pdf(paste0(plo_dir,"Biotypes_full", outname, ".pdf"), height= 4, width=6)
print(X) 
dev.off()
print("done")