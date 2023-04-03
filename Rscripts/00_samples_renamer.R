# Renames the raw_count columns to have better sample names there
# input:
#  - metadata : column 'former_name' has the ancient names, 
#               column  'sample' has the new names we want for the samples
#  - raw_counts: the columns perfectly matching with 'former_name', 
# output:
#  - raw_counts with the names of columns perfectly matching 'sample' metadata
#
# usage:
# Rscript tflex/Rscripts/00_samples_renamer.R config_postprocess.yml


# @author : johaGL 2023

library(tidyverse)
library(yaml)

############### variables 

args = commandArgs(trailingOnly=TRUE)
configfn <- args[1] # config file

# reading config parameters : 
config <- yaml::read_yaml(configfn)
mywdir <- config$mywdir
metadatafile <- config$metadata
coufi <- config$counts

setwd(mywdir)

cou_df = read.table(coufi, header=TRUE, sep='\t', check.names=FALSE)
metadf = read.table(metadatafile, header=TRUE, sep=',')

cou_df = tibble::column_to_rownames(cou_df, "ensembl_id")


intersecti_sa = intersect(colnames(cou_df), metadf$sample) 
if (length(intersecti_sa) == length(colnames(cou_df))){
  cat("\n")
  print("no need to rename the raw counts columns, already match metadata$sample")
}else{
    # construct a key-value object (dictionnary, in R is list)
    kvobj = list()
    
    for (i in 1:nrow(metadf)){
      theformer = metadf$former_name[i]
      thenew = metadf$sample[i]
      kvobj[[theformer]] <- thenew
    }
    
    # build new column names
    newnames = c()
    for (nc in colnames(cou_df)){
      new_n = kvobj[[nc]]
      newnames = c(newnames, new_n)
    }

    colnames(cou_df) = newnames
    
    cou_df = tibble::rownames_to_column(cou_df, var="ensembl_id")
    write.table(cou_df, paste0(mywdir, coufi),
               sep='\t', row.names = FALSE  )
    
    print(paste("saved ", paste0(mywdir, coufi), "with new column names"))
}


