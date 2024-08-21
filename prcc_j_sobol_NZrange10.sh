#!/bin/bash
#SBATCH --job-name=prcc_j
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --output=prcc_j_%A_%a.out
#SBATCH --error=prcc_j_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anatolii.sorokin@oist.jp
#SBATCH --array=1-380 #1-5#

echo "cores:" $OMP_NUM_THREADS

cwdir=`pwd`
module load R

fname=$1
odir=$2
target=$3
N=$4
n=${SLURM_ARRAY_TASK_ID}
mkdir -p "$odir_$N"

echo "Fname: $fname"
echo "Output directory:"
echo "$odir"

wd=$PWD
#tempdir=$(mktemp -d /flash/GoryaninU/anatoly/prcc_j.XXXXXX)
#echo $tempdir

# enter the temporary directory
#cd $tempdir
#cp /bucket/GoryaninU/FBA_sensitivity/prcc_*.R ./

Rscript prcc_j_sobol_NZrange.R $fname $odir $target 10 $n $N

sync
sleep 3

