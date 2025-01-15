#!/bin/bash
# set the partition where the job will run
#SBATCH --partition=bigmem

# set the number of cores
#SBATCH --cpus-per-task=96
# ask for more memory than 2Gb
#SBATCH --mem=500G

# Do one after the other

module load minimap2/2.24
source /soft/system/software/Miniconda3/20230719/bin/activate nanopolish

# S. uvarum

BARCODES=("FAO88471" "FAO83464" "FAO86487" "FAO79001")
DUMP="Outputs/EvolutionNanopore/nanopolish/Suva"
mkdir -p Inputs/ALBAMAR_01/nanopolish
mkdir -p Outputs/EvolutionNanopore/nanopolish/polyA_tail
mkdir -p $DUMP

for SAMPLE in ${BARCODES[@]}; do
    FAST5DIR=$DUMP"/"$SAMPLE"_1.fast5"
    FQ="/datasets/dRNA/Nanopore_yeasts_12_2020/From_HDD/ALBAMAR_01/"*$SAMPLE"/"$SAMPLE"_1_0.fastq.gz"
    ZFAST5="/datasets/dRNA/Nanopore_yeasts_12_2020/From_HDD/ALBAMAR_01/"*$SAMPLE"/"$SAMPLE"_1_0.fastq.gz"
    SUMMARY="/datasets/dRNA/Nanopore_yeasts_12_2020/From_HDD/ALBAMAR_01/"*$SAMPLE"/sequencing_summary.0.txt"
    GENOME="/datasets/CMonta/EvolutionNanopore/Inputs/Genomes/Suva_genome.fna"
    SAMPREFIX=$DUMP"/"$SAMPLE

    mkdir -p $DUMP"/"$SAMPLE
    # Write the bash script

    echo "#!/bin/bash" > $SAMPLE"_polya.sh"
    echo "#SBATCH --partition=bigmem" >> $SAMPLE"_polya.sh"
    echo "#SBATCH --cpus-per-task=20" >> $SAMPLE"_polya.sh"
    echo "#SBATCH --mem=90G" >> $SAMPLE"_polya.sh"
    echo "#SBATCH --job-name=${SAMPLE}" >> $SAMPLE"_polya.sh"
    echo "#SBATCH -o ./'slurm-'%x'-'%j'.out'" >> $SAMPLE"_polya.sh"
    echo "" >> $SAMPLE"_polya.sh"

    echo "FAST5DIR="$FAST5DIR >> $SAMPLE"_polya.sh"
    echo "FQ="$FQ >> $SAMPLE"_polya.sh"
    echo "SUMMARY="$SUMMARY >> $SAMPLE"_polya.sh"
    echo "GENOME="$GENOME >> $SAMPLE"_polya.sh"
    echo "SAMPREFIX="$SAMPREFIX >> $SAMPLE"_polya.sh"
    echo "" >> $SAMPLE"_polya.sh"

    echo "module load minimap2/2.24" >> $SAMPLE"_polya.sh"
    echo "module load SAMtools" >> $SAMPLE"_polya.sh"
    echo "source /soft/system/software/Miniconda3/20230719/bin/activate nanopolish" >> $SAMPLE"_polya.sh"
    echo "" >> $SAMPLE"_polya.sh"

# Extract compressed files
    echo "echo 'extract compressed fast5'" >> $SAMPLE"_polya.sh"
    echo "tar -xf "$(echo "/datasets/dRNA/Nanopore_yeasts_12_2020/ALBAMAR_01/"*$SAMPLE"/"$SAMPLE"_1.fast5.tar.gz")" -C "$DUMP"/"$SAMPLE"/" >> $SAMPLE"_polya.sh"
    echo "mv "$DUMP"/"$SAMPLE"/fast5 "$FAST5DIR >> $SAMPLE"_polya.sh"

# Index fastq

    echo "echo 'indexing fastq'" >> $SAMPLE"_polya.sh"
    echo "nanopolish index -d $FAST5DIR "$(echo $FQ)" -s "$(echo $SUMMARY)" --verbose" >> $SAMPLE"_polya.sh"

# Mapping part

    echo "echo 'Mapping'" >> $SAMPLE"_polya.sh"
    echo "minimap2 -t 20 -ax splice -uf -L --eqx -k14 --secondary=no -G 3000 $GENOME $FQ > ${SAMPREFIX}.sam" >> $SAMPLE"_polya.sh"
    echo "samtools view -S -b ${SAMPREFIX}.sam > ${SAMPREFIX}.bam" >> $SAMPLE"_polya.sh"
    echo "samtools sort ${SAMPREFIX}.bam -o ${SAMPREFIX}.sorted.bam" >> $SAMPLE"_polya.sh"
    echo "samtools index ${SAMPREFIX}.sorted.bam" >> $SAMPLE"_polya.sh"
    echo "rm -rf ${SAMPREFIX}.bam" >> $SAMPLE"_polya.sh"

# Run polya pipeline
    echo "echo 'Estimating polyA tail'" >> $SAMPLE"_polya.sh"
    echo "nanopolish polya --threads=20 --reads="$(echo $FQ)" --bam=${SAMPREFIX}.sorted.bam --genome=$GENOME > Outputs/EvolutionNanopore/nanopolish/${SAMPLE}_polya_results.tsv" >> $SAMPLE"_polya.sh"

    sbatch $SAMPLE"_polya.sh"

done
