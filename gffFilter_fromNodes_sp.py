import sys
import pandas as pd

try:
    functionName, GFFfile, SummaryFile, Node, output = sys.argv
except:
    print("error: Calculate_transition_matrix.py <GFF file path> <Path to summary file> <Desired Node to study> <Output gff file>")
    #quit()

#GFFfile = "/home/jmontanes/Documents/0-Important_files/Yeasts/VCF_files/Parsing_SNPs/ORF_GFF/s_cerevisiae_annotations_CDS_v64_1_1_liftOverdone.gff"
#SummaryFile = "/home/jmontanes/Documents/IQtree_Gene_duplication/Yeasts/Standard_annotation/OutputsR/Summary/YEASTS_gene_summaryshort.tsv"
#output = "Potato.txt"

Summary_table = pd.read_table(SummaryFile, sep='\t', comment="#", dtype='str')
Transcripts = Summary_table[Summary_table["AgeNode"] == Node].Transcript
outputFile = open(output, "w+")
with open(GFFfile, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        if attr.startswith("Parent="):
            transcript = attr.split("=")[1].split("\n")[0]
            if sum(Transcripts == transcript) > 0:
                outputFile.writelines(line)

        #if gxfToExtract_path.split(".")[-1] == "gtf":

        #elif gxfToExtract_path.split(".")[-1] == "gff" or gxfToExtract_path.split(".")[-1] == "gff3":
outputFile.close()

