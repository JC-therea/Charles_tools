#!/bin/bash
# set the partition where the job will run
#SBATCH --partition=normal,bigmem,long

# set the number of cores
#SBATCH --cpus-per-task=12
# ask for more memory than 2Gb
#SBATCH --mem=20G

#SBATCH --job-name=featureCounts
#SBATCH -o ./"slurm-"%x"-"%j".out"

echo "To run the program"
echo "wd=/path/to/mapped/files"
echo "FILENUMBER=$(ls ${wd} 2>/dev/null | wc -l )"
echo "sbatch $wd --array=1-$FILENUMBER src/fcOnShiva.sh"

echo $2
