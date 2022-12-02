
import sys

try:
    functionName, filePath, featureChosen, outFilePath = sys.argv
except:
    print("error: ribORFgenePred2gtf.py <ribORF genePred file> <ribORF features splitted by comma or All> <Output file>")
    quit()

featureChosenList = featureChosen.split(",")
for feature in featureChosenList:
    if feature not in ["All","canonical","extension","external","internal","noncoding","overlap.uORF","polycistronic","readthrough","truncation","uORF"]:
        print("error: Not correct feature type " + feature)
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
       ORFexon = 1
       attribute_out = 'gene_id "' + geneID.split("-T")[0] + '"; transcript_id "' + geneID + '"; ID "' + ORF + "-" + str(ORFexon) + '"; Parent "' + geneID + '";'
       NewORFstart = int(ORFstart) + 1
       if ORFType in featureChosenList or "All" in featureChosenList:
           if int(exonNumber) == 1:
               CDS_text = chr + '\t' + "ribORF" + '\t' + "ORF" + '\t' + str(NewORFstart) + '\t' + ORFend + '\t' + '.' + '\t' + ribORFstrand + '\t' + "0" + '\t' + attribute_out + "\n"
               outFile.writelines(CDS_text)

           else:
               CDS_text = chr + '\t' + "ribORF" + '\t' + "ORF" + '\t' + str(NewORFstart) + '\t' + ORFend + '\t' + '.' + '\t' + ribORFstrand + '\t' + "0" + '\t' + attribute_out + "\n"
               ORF = geneID + "_ORFID=" + nORF
               exonStarts_List = exonStarts.split(",")
               exonEnds_List = exonEnds.split(",")
               ORFStart_TrStarts = 0
               ORFEnd_TrStarts = 0
               # Locate the start and the end of the ORF
               for exonStarts_element in exonStarts_List[0:(len(exonStarts_List) - 1)]:
                   if int(exonStarts_element) < NewORFstart:
                       ORFStart_TrStarts += 1
                   if int(exonStarts_element) < int(ORFend):
                       ORFEnd_TrStarts += 1
               lines_to_print = ORFStart_TrStarts - ORFEnd_TrStarts

               if lines_to_print == 0:

                   if ORFType in featureChosenList:
                       outFile.writelines(CDS_text)
               else:
                   for exon in range(ORFStart_TrStarts - 1, ORFEnd_TrStarts):
                       attribute_out = 'gene_id "' + geneID.split("-T")[0] + '"; transcript_id "' + geneID + '"; ID "' + ORF + "-" + str(ORFexon) + '"; Parent "' + geneID + '";'
                       ORFexon += 1
                       if exon == (ORFStart_TrStarts - 1):
                           CDS_text = chr + '\t' + "ribORF" + '\t' + "ORF" + '\t' + str(NewORFstart) + '\t' + exonEnds_List[exon] + '\t' + '.' + '\t' + ribORFstrand + '\t' + "0" + '\t' + attribute_out + "\n"
                           outFile.writelines(CDS_text)
                       elif exon == (ORFEnd_TrStarts - 1):
                           CDS_text = chr + '\t' + "ribORF" + '\t' + "ORF" + '\t' + exonStarts_List[exon] + '\t' + ORFend + '\t' + '.' + '\t' + ribORFstrand + '\t' + "0" + '\t' + attribute_out + "\n"
                           outFile.writelines(CDS_text)
                       else:
                           CDS_text = chr + '\t' + "ribORF" + '\t' + "ORF" + '\t' + exonStarts_List[exon] + '\t' + exonEnds_List[exon] + '\t' + '.' + '\t' + ribORFstrand + '\t' + "0" + '\t' + attribute_out + "\n"
                           outFile.writelines(CDS_text)

outFile.close()
