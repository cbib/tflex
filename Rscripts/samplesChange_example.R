# johaGL 2022
library(tidyverse)
## NOTE: former matrices and metadata are in zip files
mydir <- "~/cocultureProj/"
setwd(mydir)
# ------------ open all metadata and raw counts files
# suffixes in metadata and raw counts
suffixes <- c("Human","Rat")
mydfs <- list()
mycounts <- list()
openmetadata <- function(ubi){
  oh <- read.table(paste0( ubi,".csv"),sep=",", header=TRUE)
  return(oh)
}
opencounts <- function(ubi){
  oh <- read.table(paste0( ubi,".tsv"),sep="\t", header=TRUE,
                   check.names = FALSE)
  return(oh)
}
for(suf in suffixes){
  mydfs[[suf]] <- openmetadata(paste0("metadata/samples",suf))
  mycounts[[suf]] <- opencounts(paste0("datarc/raw_counts",suf))
}

# ------------change names in rat----------------
mydfs$Rat$sample
mydfs$Rat$formername <- mydfs$Rat$sample
tmp <- unname(sapply(mydfs$Rat$sample,function(x){
  str_replace(x,"Neurones-avec-insert", "Neu-insert")}))
tmp <- unname(sapply(tmp,function(x){
  str_replace(x,"Neurones", "Neu-control-insert")}))
tmp <- unname(sapply(tmp,function(x) unlist(str_split(x, "_"))[2]))
# I encounter excel auto-column problem P3, P4 P5 ....
# > tmp
# [1] "Neu-control-insert" "Neu-control-insert" "Neu-control-insert"
# [4] "Neu-control-insert" "Neu-insert-P3"      "Neu-insert-P4"     
# [7] "Neu-insert-P5"      "Neu-insert-P6"   
tmp2 <- c()
for (i in tmp){
  if (str_starts(i,"Neu-insert-") ){
    # str_replace(i, "[3456]", "")
    tmp2 <- c(tmp2,"Neu-insert")
  }else{tmp2 <- c(tmp2,i)}
}
reps <- rep(c('1','2','3','4'),2)
oknamesrat <- paste0(tmp2,"-",reps)
mydfs$Rat$sample <- oknamesrat
View(mydfs$Rat) # ok everything, 
# TODO: apply to count matrix
colnames(mycounts$Rat)
# DO NOT TOUCH first colname ! 
if (all(colnames(mycounts$Rat)[2:length(colnames(mycounts$Rat))] ==
        mydfs$Rat$formername)){
  print("ok perfect match, assigning new colname to matrix")
  colnames(mycounts$Rat)[2:length(colnames(mycounts$Rat))] <- mydfs$Rat$sample
}else {print("ERROR!!! impossible to change colnames with this easy method")}
# --------------------change names in human-------------------
mydfs$Human$sample
mydfs$Human$formername <- mydfs$Human$sample
tmpH <- mydfs$Human$sample
tmpH <- unname(sapply(tmpH,function(x){
  str_replace(x,"P3-cultives-en-NBM_glutamine", "P3-nbmGln")}))
tmpH <- unname(sapply(tmpH,function(x){
  str_replace(x,"P3-cultives-en-NBM", "P3-nbm")}))
tmpH <- unname(sapply(tmpH,function(x){
  str_replace(x,"P3-en-insert-sur-Neurones", "P3-insert")}))
tmpH <- unname(sapply(tmpH,function(x) unlist(str_split(x, "_"))[2]))
tmpH[1:3] <- rep("P3-control-insert",3)
repsH <- c(c('1','2','3'), rep(c('1','2','3','4'),3))

mydfs$Human$sample <- paste0(tmpH,"-",repsH)
# TODO : apply to count matrix
if (all(colnames(mycounts$Human)[2:length(colnames(mycounts$Human))] ==
        mydfs$Human$formername)){
  print("ok perfect match, putting new names in matrix")
  colnames(mycounts$Human)[2:length(colnames(mycounts$Human))] <- mydfs$Human$sample
}else {print("ERROR!!! impossible to change colnames with this easy method")}

# ------------ save all metadata and counts with new names
savemetadata <- function(thing,ubi){
  write.table(thing,paste0( ubi,"-new.csv"),sep=",", col.names=TRUE, row.names=FALSE)
}
savecounts <- function(thing,ubi){
  write.table(thing,paste0( ubi,"-new.tsv"),sep="\t", col.names=TRUE, row.names = FALSE)
}
for (suf in suffixes){
  savemetadata(mydfs[[suf]], paste0("metadata/samples",suf))
  savecounts(mycounts[[suf]], paste0("datarc/raw_counts",suf))
}

