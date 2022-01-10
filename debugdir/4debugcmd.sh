#!/bin/bash
########################## Slurm options #####################
#SBATCH --job-name=debugIII
#SBATCH --workdir=/home/jgalvis/tflex/debugdir/
#SBATCH --mail-user=juana7@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --exclusive

################################################################

scontrol show job $SLURM_JOB_ID

echo "testing slurm alone"
sleep 8
touch debug_simpler.txt
echo "bye bye"

