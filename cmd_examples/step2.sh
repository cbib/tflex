#!/bin/bash
########################## Slurm options
#SBATCH --job-name=snakestar
#SBATCH --output=/mnt/cbib/thesis_gbm/mubriti_202303/scr1/slurm_output/snakestar_%j.out
#SBATCH --workdir=/mnt/cbib/thesis_gbm/mubriti_202303/scr1/
#SBATCH --mail-user=juana7@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00
#SBATCH --exclusive
##################################################

scontrol show job $SLURM_JOB_ID

tflexPath="/mnt/cbib/thesis_gbm/tflex"

configpath="/mnt/cbib/thesis_gbm/mubriti_202303/scr1/config_mapping.yml"

module load snakemake
module load fastp
module load multiQC
module load STAR/2.7.10a # new Star

# 2 #
# Calling mapping with snakemake, separately for both species
# note: --j for threads ok, --cores 1 to force sequentially run over samples
echo "running mapping "
snakemake -s $tflexPath/Snake_starmp.smk --cores 1 -j 12 \
    --configfile $configpath --latency-wait=30


