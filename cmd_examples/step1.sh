#!/bin/bash
########################## Slurm options
#SBATCH --job-name=snakefp
#SBATCH --output=/mnt/cbib/thesis_gbm/mubriti_202303/scr1/slurm_output/snakefp_%j.out
#SBATCH --workdir=/mnt/cbib/thesis_gbm/mubriti_202303/scr1/
#SBATCH --mail-user=juana7@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --ntasks-per-node=12
#SBATCH --time=12:00:00
#SBATCH --exclusive
##################################################

scontrol show job $SLURM_JOB_ID

tflexPath="/mnt/cbib/thesis_gbm/tflex"
configpath="/mnt/cbib/thesis_gbm/mubriti_202303/scr1/config_mapping.yml"

module load snakemake
module load fastp
module load multiQC


# 1 # 
# Calling Snake_fastp.smk : qc with multiqc, trimming + adapter removal by fastp

# avoid message "locked directory" : --unlock option
snakemake --unlock -s $tflexPath/Snake_fastp.smk --cores 1 \
 --configfile $configpath

# run
snakemake -s $tflexPath/Snake_fastp.smk --cores 12 \
 --configfile $configpath --latency-wait=10


# 2 #
