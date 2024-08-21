#!/bin/bash
#SBATCH --job-name=cFlRanges
#SBATCH --partition=compute
#SBATCH --time=74:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --output=cFlRanges_%j.out
#SBATCH --error=cFlRanges_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anatolii.sorokin@oist.jp

echo "cores:" $OMP_NUM_THREADS

cwdir=`pwd`
module load R

mdir=$1
odir=$2
ndir=$3
fnum=$4

echo "arguments:" $mdir $ndir $fnum

Rscript checkFluxRanges.R $mdir $odir $ndir $fnum
