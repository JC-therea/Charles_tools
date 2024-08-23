#!/bin/bash
## set the partition where the job will run
#SBATCH --partition=bigmem

# set the number of cores
#SBATCH --cpus-per-task=16
# ask for more memory than 2Gb
#SBATCH --mem=90G

# Job name
#SBATCH --job-name=FunannotateUpdate-BBQ

#SBATCH -o ./"slurm-"%x"-"%j".out"

module load Miniconda3/4.9.2
source /soft/system/software/Miniconda3/4.9.2/bin/activate funannotate
module load SAMtools/1.12-GCC-10.2.0
module load Perl/5.32.0-GCCcore-10.2.0

export FUNANNOTATE_DB=/datasets/CMonta/funannotate_db/

Y="BBQ"

GENOME="/datasets/dRNA/Gauthier/24052024/genomes/"$Y".Final.fasta"
GENOMEFIX="Output/Funannotate/"$Y"/"$Y_"genome.sorted.fna"
GFF="/projects_eg/projects/jmontanes/Colaborations/Gauthier/Outputs/customToga/"$Y".toga.mod.fun.gff"
dRNA="/projects_eg/projects/jmontanes/Colaborations/Gauthier/Outputs/Funannotate/"$Y"/rnabloom/rnabloom.transcripts.fa"

mkdir -p "Output/Funannotate/"$Y"/"

# Fix genome 

#sed 's/ .*//g' $GENOME > $GENOMEFIX

funannotate update -f $GENOMEFIX -g $GFF -o $Y"Bloom" \
--nanopore_mrna $dRNA --no_trimmomatic --pasa_db mysql \
--stranded F --jaccard_clip --species "Saccharomyces cerevisiae" \
--cpus 16 --no-progress --max_intronlen 1200
