#!/bin/bash
#SBATCH --job-name=cSobEcoFBAd
#SBATCH --partition=compute
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --output=cSobEcoFBA_%A_%a.out
#SBATCH --error=cSobEcoFBA_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anatolii.sorokin@oist.jp
#SBATCH --array=1-16 #1-64 #1-512 #1-5#

echo "cores:" $OMP_NUM_THREADS
echo "array:" ${SLURM_ARRAY_TASK_ID}

cwdir=`pwd`

chunkSize=8192
model=$1
odir=$2
seed=255

mkdir -p odir 
export model odir chunkSize seed
echo $SLURM_ARRAY_TASK_ID $model $odir $chunkSize $seed

module load python/3.10.2

python3 calc_FBA_ecoSobol_dump.py

sync
sleep 3

