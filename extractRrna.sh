genome=$1

mkdir -p jmontanes

source /soft/system/software/Miniconda3/20240927/bin/activate barrnap
barrnap --kingdom euk --threads 18 --outseq "jmontanes/rRNA.barrnap.fa" $genome > jmontanes/rRNA.gff
