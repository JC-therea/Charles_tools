
import sys

try:
    functionName, filePath, outFilePath = sys.argv
except:
    print("error: ribORFgenePred2gtf.py <ribORF genePred file> <Output file>")
    quit()


ORF_info_path = filePath
outPath = outFilePath + ".gtf"

# Get coordinates in sequence
outFile = open(outPath, "w+")
with open(ORF_info_path, "r") as file:
   for line in file:
       ribORF, chr, strand, transcriptStart, TranscriptEnd, ORFstart, ORFend, exonNumber, exonStarts, exonEnds = line.split("\t")
       geneID, ribORFchr, strand_nORF_transcriptLength, ribORFstart, ribORFend_ORFType_StartCodon = ribORF.split(":")
       ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
       ribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")
       ORF = geneID + "_ORFID=" + nORF

       attribute_out = 'gene_id "' + geneID.split("-T")[0] + '"; transcript_id "' + geneID + '"; ID "' + ORF + '"; Parent "' + geneID + '";'

       NewORFstart = int(ORFstart) + 1
       CDS_text = chr + '\t' + "ribORF" + '\t' + "ORF" + '\t' + str(NewORFstart) + '\t' + ORFend + '\t' + '.' + '\t' + ribORFstrand + '\t' + "0" + '\t' + attribute_out + "\n"
       outFile.writelines(CDS_text)

outFile.close()
