#!/usr/bin/env Rscript
# johaGL 2022
library(tidyverse)
library(msigdbr)
library(fgsea)
library(data.table)

#save n top pathways: 
ntops = 20

args = commandArgs(trailingOnly=TRUE)
#gv <- yaml::read_yaml("~/cocultureProj/scr2/set1_Neu_insert-vs-control.yaml")
gv <- yaml::read_yaml(args[1])

reqcontrast <- gv$requiredcontrast 
outname <- gv$outname
equi_id <- gv$equivalentid  # inspect gene id's  ('.':ensembl_gene_id_version) !
shortname <- gv$shortname
trustedcolumn <- gv$trustedcolumn
#SPECIESensmbldsetname <- gv$SPECIESensmbldsetname
# conditionLEVELS <- gv$conditionLEVELS
chosencollections <- c("GO:BP", "CP:KEGG")
print(paste("internal list fullspeciesname, key invoked by 'shortname' ! :",
            "MANDATORY 'mouse' OR 'rat' OR 'human'"))
fullspeciesname <- list("mouse"="Mus musculus", "rat"="Rattus norvegicus",
                    "human"="Homo sapiens")

source(paste0(args[2], "Rscripts/func.R"))

# START
setwd(gv$mywdir)
o_dir <- paste0(gv$mywdir,"results_",outname,"/")
resdirs <- getoutputdirs(outer=o_dir)
rds_dir = resdirs[1]; tabl_dir = resdirs[2]; plo_dir = resdirs[3] 
strcontrast = paste(reqcontrast[2:3],collapse="_vs_")

degs <- read.table(paste0(tabl_dir,"DEtop_",strcontrast,"_",outname,".tsv"),
header=TRUE,  sep='\t')

# -- smal pretreatment, compulsory
print(paste("saving true list sent to GSEA, as symbols must not be repeated",
            ", nor contain dot '.' (used when make.unique in 1_prepdata)"))

truelist <- degs[,c("ensembl_id", "pvalue","log2FoldChange", trustedcolumn)]
awaydots <- sapply(truelist[[trustedcolumn]], function(x){
  return(unname(unlist(str_split(x,"\\."))[1]))
})  
truelist[[trustedcolumn]] <- awaydots # no longer unique, pick min pval if multiplets
truelist.true <- truelist %>% group_by_at(trustedcolumn) %>% 
  slice_min(pvalue, with_ties = FALSE)
rm(truelist)
write.csv(truelist.true, 
        paste0(tabl_dir,"fgseaquery_",strcontrast,"_",outname,".csv"), 
        row.names=FALSE) # csv saves header by default

minnonzero <- min(truelist.true$pvalue[truelist.true$pvalue > 0])
truelist.true$pvalue[truelist.true$pvalue==0] <- minnonzero

truelist.true <- truelist.true %>% 
	mutate(mystat = -1*(log10(pvalue))*(log2FoldChange/abs(log2FoldChange))
	)
veclg <- truelist.true$mystat
names(veclg) <- truelist.true$symbol_unique

head(veclg)
tail(veclg)
# -- end pretreatment

for (kk in chosencollections){
  print(paste("*    using collection:",kk))
  heregmt <- givegmt(kk,fullspeciesname[[shortname]])
  msigdbr_list =split(x=heregmt$gene_symbol, f=heregmt$gs_name)
  set.seed(42)
  fgseaRes <- fgsea(pathways=msigdbr_list,
                    stats=veclg,
                    minSize=3,
                    maxSize=Inf, 
                    scoreType="std",
                    nperm=100000)
  saveRDS(fgseaRes, paste0(rds_dir,"fgsea-",kk,"_",strcontrast,"_",outname,".rds"))
  fullfn <- paste0(tabl_dir, "fgseafull-",kk,"_",strcontrast,"_",outname,".tsv")
  fwrite(fgseaRes, file=fullfn, sep="\t", sep2=c("", " ", ""))
  tata <- read.table(fullfn, header=TRUE, sep="\t")
  upo <- tata %>% filter(NES > 0) %>% arrange(padj) 
  dow <- tata %>% filter(NES < 0 ) %>% arrange(padj)
  savetsv(head(upo,n=ntops), paste0(tabl_dir, "fgseaUP-",kk,"_",strcontrast,"_",outname,".tsv"))
  savetsv(head(dow, n=ntops),paste0(tabl_dir, "fgseaDOWN-",kk,"_",strcontrast,"_",outname,".tsv"))
  #fgseaRes <- readRDS(paste0(rds_dir,"fgsea-",kk,"_",strcontrast,"_",outname,".rds"))
  
  head(fgseaRes[order(pval),])
  
  paplo <- paths4plot(fgseaRes, 10)
  pdffin <- paste0(plo_dir,"fgsea-",kk,"_",strcontrast,"_",outname,".pdf") 
  pdf(pdffin, width=13, height=5)
  print(plotGseaTable(msigdbr_list[c(paplo$up,rev(paplo$down))],
                veclg, fgseaRes))
  dev.off()
}
