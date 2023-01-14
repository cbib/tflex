#!/usr/bin/env Rscript
# example: Rscript --vanilla 0_verify_options.R ratvars.yaml
args = commandArgs(trailingOnly=TRUE)
gv <- yaml::read_yaml(args[1])
setwd(gv$mywdir)

metadatafn <- gv$metadatafile
cfn <- gv$countsfile
conditionLEVELS <- gv$conditionLEVELS
reqcontrast <- gv$requiredcontrast 
outname <- gv$outname
equi_id <- gv$equivalentid
shortname <- gv$shortname
samplestodrop <- gv$samplestodrop
print(metadatafn)
print(cfn)
print(conditionLEVELS)
print(reqcontrast)
print(outname)
print(equi_id)
print(shortname)
print(samplestodrop)
