# Charles_tools
Here I will keep all those standalone scripts that are useful for bioinformatics and for me.  
In the directory singleUse are those scripts that were useful for a specific purpose but probably will not be very useful in general.

Here goes some explanation about the others:

## AnnotationQC.py

This script works like this:   
`AnnotationQC.py <Input file> <Output file>`  
From a gff/gtf file, this script extracts the sizes per transcript of the CDSs, UTRs, introns, etc. Output file is in tsv format.

## extract_nucleotides.py

From a fasta file and an annotation file this script extract the nucleotide sequence. Work like this:  
`extract_nucleotides.py -gff $gff -fasta $genome -f <features selected (3rd column)> -o <output file in fasta format>`

## gffCompareExtraction.py

Script to extract specific gffCompare identifiers from a gff file. It works like this:
`gffCompareExtraction.py <gffCompare.tracking file> <gffCompare identifiers chosen separated by ','> <gff to extract sequences> <Output gff>`

## InSilico_translation.py

Get all the possible ORFs from a multifasta file. The minimum size of the ORFs has to be indicated when running the script.   
`inSilico_translation.py <Directory with multifasta files> <Amino acid threshold>`

## msaOrfs.py

This program uses genepred, orthology and msa files to detect orfs that are aligned. This program will go to singleUse directory when I finish my current manuscript. 
`msaOrfs.py -i <genePref Files splitted by commas> -orth <orthology table from proteinOrtho> -msa <directory where all the MSA were done by transcript> -o <output file>`

## NewickTreeView.py

Simple tool to visualize newick trees.
`NewickTreeView.py newick_tree_file`

## PatternReporter.py

Report the number of times a pattern appears in a given multifasta file per transcript.  
`PatternReporter.py -p <NUcleotide pattern to search> -fasta <Path to the multifasta file> -o <Output file in fasta format>`

## singletReport.py and TripletReport.py

Programs to report the nucleotide composition per fasta in multifasta files. singletReport reports each single nucleotide and TripletReport does the same but for triplets.
`singletReport.py -fasta <Path to the multifasta file> -o <Output file>`
`TripletReport.py -fasta <Path to the multifasta file> -o <Output file>`

## transcriptExtender.sh

Simple pipeline to extend transcript from a gff file. It requires agat.
`transcriptExtender.sh <input.gff> <out directory>`

