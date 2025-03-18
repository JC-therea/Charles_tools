inputDir=$1

files=$(ls ${inputDir}/*.fastq ${inputDir}/*.fq ${inputDir}/*.fastq.gz ${inputDir}/*.fq.gz 2>/dev/null)

for fq in ${files}; do
	reads=$(echo $(zcat $fq |wc -l)/4|bc)

	echo -e $(basename "$fq")'\t'$reads
done
