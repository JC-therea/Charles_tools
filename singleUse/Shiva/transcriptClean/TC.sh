#!/bin/bash
# set the partition where the job will run
#SBATCH --partition=normal,long,bigmem

# set the number of cores
#SBATCH --cpus-per-task=24
# ask for more memory than 2Gb
#SBATCH --mem=90G

#SBATCH --job-name=TranscriptClean
#SBATCH -o ./"slurm-"%x"-"%j".out"

# To run the second cleaning method with TranscriptClean
module load minimap2
module load Miniconda3/4.9.2
module load SAMtools

### Path configuration

FQDIR="/projects_eg/projects/jmontanes/Colaborations/CPapadopoulos/Outputs/fmlrc/"
TRANSCRIPTCLEAN_REPO="TranscriptClean"

##########################
##########################

########## File

FASTQFILES=($(ls -1 ${FQDIR}/*.fasta 2>/dev/null ))
i=$(($SLURM_ARRAY_TASK_ID - 1)) ## bash arrays are 0-based
THISFASTQFILE=${FASTQFILES[i]}
ID=$(echo $(basename "$THISFASTQFILE") | sed -E 's/_fmlrc\.RNA\.fasta$//')

mkdir -p "Outputs/TranscriptClean/"$ID"/"

# Select the right genome

if [ "$ID" == "Spar" ]; then
    echo "Spar"
    GENOME="/projects_eg/projects/chris/Enriched_strains_RNAseq/Strains_Genomes_Comparisons/Genomes_T2T_fasta/PAR/PAR/Saccharomyces_paradoxus_CBS432.fasta"

elif [ "$ID" == "Suva" ]; then
    echo "ID is Suva"
    GENOME="/projects_eg/projects/chris/Enriched_strains_RNAseq/Strains_Genomes_Comparisons/Genomes_T2T_fasta/UVA/CBS/Saccharomyces_uvarum_CBS7001.fna"

elif [ "$ID" == "Scer" ]; then
    echo "ID is Scer"
    GENOME="/projects_eg/projects/chris/Enriched_strains_RNAseq/Strains_Genomes_Comparisons/Genomes_T2T_fasta/CER/R64/saccharomyces_cerevisiae_R64-1-1.fsa"
fi

SAMPREFIX_fmlrc="Outputs/"$ID
SAMPREFIX_TC_WD="Outputs/TranscriptClean/"$ID"/"
SAMPREFIX_TC="Outputs/TranscriptClean/"$ID"/"$ID

mkdir -p $SAMPREFIX_TC_WD

TCBAM=$SAMPREFIX_TC".sorted.bam"
OUTBAM=$SAMPREFIX_TC".MAPQ.sorted.bam"
OUTFASTA=$SAMPREFIX_TC".clean.fa"

##########################
##########################

############################################
############## Main script #################
############################################

echo "Processing "$ID

if [ -f "$OUTFASTA" ]; then
    echo "File $OUTFASTA exist. Stopping the script."
    exit 1
fi

# Mapping
echo "Starting mapping"
echo "Genome "$GENOME
echo "file "$THISFASTQFILE

minimap2 -t 24 -ax splice -uf -L -k14 --secondary=no -G 2600 $GENOME $THISFASTQFILE > $SAMPREFIX_fmlrc".sam"
samtools sort $SAMPREFIX_fmlrc".sam" -o $SAMPREFIX_fmlrc".sorted.sam"

# Cleaning with transcript clean

module load Miniconda3/20230719

echo "Starting transcriptClean"

python $TRANSCRIPTCLEAN_REPO"/"TranscriptClean.py --sam $SAMPREFIX_fmlrc".sorted.sam" --genome $GENOME --outprefix $SAMPREFIX_TC_WD -t 24
mv $SAMPREFIX_TC_WD"/TC_clean.sam" $SAMPREFIX_TC".sam"

##### Sort BAM file #####

module load SAMtools/1.12-GCC-10.2.0

echo "Sort and index the bam file"

samtools view -S -b $SAMPREFIX_TC".sam" > $SAMPREFIX_TC".bam"
samtools sort $SAMPREFIX_TC".bam" -o $SAMPREFIX_TC".sorted.bam"
samtools index $SAMPREFIX_TC".sorted.bam"
rm -rf $SAMPREFIX_TC".bam"

##### MAPQ filter #####
echo "MAPQ filter of the bam file"

samtools view -Sb -q 5 -f 0,16 -h $TCBAM  > $OUTBAM
samtools fasta -f 0,16 $OUTBAM > $OUTFASTA
