#!/bin/sh

# Here is a try cath statement
CONF="genomicCoverageConfig.txt"
if [ -f "$CONF" ];
then
	source $CONF
	echo "Config file found"
else
	echo "Config file not detected, making..."
	echo "RNABAM=" > genomicCoverageConfig.txt
	echo "RIBOBAM=" >> genomicCoverageConfig.txt
	echo "GTF=" >> genomicCoverageConfig.txt
	echo "SparKsource=" >> genomicCoverageConfig.txt
	echo "Done! Please fill the configuration file :D"
	exit 0
fi

while getopts g:o: flag
do
    case "${flag}" in
        g) GENE=${OPTARG};;
        o) OUTPUT=${OPTARG};;
    esac
done

echo $GENE
echo $OUTPUT

CHR=$(awk '$3 == "gene"' $GTF | grep ${GENE}'"' | cut -f1)
START=$(( $(awk '$3 == "gene"' $GTF | grep ${GENE}'"' | cut -f 4)  ))
END=$(( $(awk '$3 == "gene"' $GTF | grep ${GENE}'"' | cut -f 5)  ))

echo "Chromosomic region $CHR:$START-$END"

python $SparKsource -cf $RNABAM $RIBOBAM -pr $CHR":"$START"-"$END \
-o $OUTPUT -gtf $GTF \
-gl RNA Ribo
