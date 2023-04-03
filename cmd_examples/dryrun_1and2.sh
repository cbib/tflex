#!/bin/bash
# DRYRUN does not runs for real
# only shows what will be done

module load snakemake

tflexPath="/mnt/cbib/thesis_gbm/tflex"

configpath="/mnt/cbib/thesis_gbm/mubriti_202303/scr1/config_mapping.yml"
# 1 #
# qc and trimming step
echo "dryrun  qc/trimming on entire dataset "


snakemake --dryrun -s $tflexPath/Snake_fastp.smk --cores 1 \
  --configfile $configpath


# 2 #

# note: --j for threads ok, --cores 1 to force sequentially run over samples
echo "dryrun mapping "
snakemake --dryrun -s $tflexPath/Snake_starmp.smk --cores 1 -j 12 \
    --configfile $configpath --latency-wait=30


