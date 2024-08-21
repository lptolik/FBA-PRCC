#!/bin/bash
#SBATCH --job-name=makeRanges
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --output=makeRanges_%j.out
#SBATCH --error=makeRanges_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anatolii.sorokin@oist.jp

echo "cores:" $OMP_NUM_THREADS

cwdir=`pwd`
module load R


mname=$1
odir=$2
fnum=64
if [ ! -z "$3" ]
  then
    fnum=$3
    echo "Fnum supplied $fnum"
fi
echo $mname, $odir, $fnum

Rscript makeRangesN.R $mname $odir $fnum

sync
sleep 3

