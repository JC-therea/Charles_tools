#!/bin/bash
# set the partition where the job will run
#SBATCH --partition=normal

# set the number of cores
#SBATCH --cpus-per-task=2
# ask for more memory than 2Gb
#SBATCH --mem=10G

#SBATCH --job-name=TranscriptClean
#SBATCH -o ./"slurm-"%x"-"%j".out"

# TRANSCRIPTCLEAN_REPO="/home/jmontanes/Documents/Software/TranscriptClean-master"

echo "HI"

FQDIR="/projects_eg/projects/jmontanes/Colaborations/CPapadopoulos/Outputs/fmlrc"

FILENUMBER=$(ls ${FQDIR}/*.fasta | wc -l )
echo "[ $(date "+%Y-%m-%d %H:%M:%S") ] Total files to analyze: $FILENUMBER"

sbatch --array=1-$FILENUMBER --parsable TC.sh
