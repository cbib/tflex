#!/bin/bash
########################## Slurm options
#SBATCH --job-name="debug"
#SBATCH --output=$HOME/tflex/debugdir/slurm_output/debug_%j.out
#SBATCH --workdir=$HOME/tflex/debugdir/
#SBATCH --mail-user=juana7@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=500M
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
##################################################

module load snakemake
module load STAR/2.7.9a


snakemake -s Sn_test.smk --cores 1 --latency-wait=15 \
--cluster-config cslurm.yaml \
--cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes}"


