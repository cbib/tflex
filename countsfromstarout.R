#!/usr/bin/env Rscript
# This is 'countsfromstarout.R' 
# takes arguments to yield counts from multiple Star output count files
# 
# @author : johaGL

library(tidyverse)
library(yaml)

############### variables 

args = commandArgs(trailingOnly=TRUE)
configfn <- args[1] # config file
outfn <- args[2] # output file
filestarcount <- "ReadsPerGene.out.tab"  # by default

# reading config parameters : mappeddir and metadata
config <- yaml::read_yaml(configfn)
mappeddir <- config$mappeddir
metadatafile <- config$metadata
metadatadf <- read.table(metadatafile, header=T, sep=",")
print("processing these samples:")
print(metadatadf)
samplesVector <- metadatadf[[1]] # first column, colname ignored

## functions

getunstrandedcounts <- function(mappeddir, sample, filestarcount){
  # note : add "/" among sample and filecountstar 
  pathforig = paste0(mappeddir,sample, "/", filestarcount )
  dforig = read.table(pathforig, header=F, sep='\t')
  rowsnoneed = c("N_noFeature", "N_unmapped","N_multimapping","N_ambiguous")
  df1 <- dforig %>% filter(! V1 %in% rowsnoneed) %>% select(V1, V2)
  colnames(df1) <- c("ensembl_id", sample)
  rownames(df1) <- df1[["ensembl_id"]]
  return(df1) # returns only ensembl_id and unstranded counts (two cols)
}

getcountsonespecies <- function(mappeddir, sampleList, filestarcount){
  dflist <- list()
  for (sample in samplesVector){
    dflist[[sample]] <- getunstrandedcounts(mappeddir, sample, filestarcount)
  }
  joined <- tibble("ensembl_id" = dflist[[1]]$ensembl_id)
  for (i in 1:length(dflist)){
    joined <- full_join(joined, dflist[[i]], by="ensembl_id") 
  }
  if (all(colnames(joined) == c("ensembl_id", samplesVector))){
    print("full_join into unique counts matrix sucessfull! Colnames are:  ")
    cat(samplesVector,sep=', ')  # just prints samples
    print("")
  }
  return(joined)  
}

## execution
if ((str_ends(configfn, "yaml") | str_ends(configfn, "yml")) && str_ends(outfn, "tsv")){
  resu <- getcountsonespecies(mappeddir, samplesVector, filestarcount)
  write.table(resu, outfn, col.names = T, row.names = F, sep='\t')
}else{print("bad files in args, did you inversed ?")}

print("DONE Rscript : 'countsfromstarout.R'")




