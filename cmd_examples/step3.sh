#!/bin/bash

# first run :
# conda activate r_env_tflex

module load snakemake
module load python/3.9.7


dir_out_counts="/mnt/cbib/thesis_gbm/mubriti_202303/raw_counts/"

out_file=${dir_out_counts}"counts_mubrit.tsv"

config_yml="/mnt/cbib/thesis_gbm/mubriti_202303/scr1/config_mapping.yml"

tflexPath="/mnt/cbib/thesis_gbm/tflex"

if [ ! -d $out_file ];then 
	mkdir -p $dir_out_counts
fi

# info mapping
python $tflexPath/extract_mapp_inf.py $config_yml $dir_out_counts

# counts 
Rscript --vanilla $tflexPath/countsfromstarout.R $config_yml $out_file
