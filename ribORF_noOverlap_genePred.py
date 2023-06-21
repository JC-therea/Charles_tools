import sys

try:
    functionName, filePath, featureChosen, ATGmode, outFilePath = sys.argv
except:
    print("error: ribORF_noOverlap_genePred.py <Input file> <Chosen features separated by , or All> <ATG mode = T/F> <Output file>")
    quit()

featureChosenList = featureChosen.split(",")
for feature in featureChosenList:
    if feature not in ["All","canonical","extension","odORF","iORF","noncoding","ouORF","dORF","readthrough","truncation","uORF"]:
        print("error: Not correct feature type " + feature)
        quit()
#filePath = "/users/genomics/jmontanes/EvolutionaryNanopore/AdditionalSamplesRiboseq/riboNovel-Scer/Scer/Correct_format_files/candidateORF.genepred.txt"
#outFilePath = "p.txt"
#featureChosen = "canonical,uORF,dORF"
#ATGmode = "T"
featureChosenList = featureChosen.split(",")
file = open(filePath,"r")
outFile = open(outFilePath,"w+")
canonicalFrame = {}
ORFstartEnd = {}

keepedORFs = {}

for line in file:
    ribORF, chr, strand, transcriptStart, TranscriptEnd, ORFstart, ORFend, exonNumber, exonStarts, exonEnds = line.split("\t")
    geneID, ribORFchr, strand_nORF_transcriptLength, ribORFstart, ribORFend_ORFType_StartCodon = ribORF.split(":")
    ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
    ribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")
    if ORFType not in featureChosenList and "All" not in featureChosenList:
        continue
    # if the gene ID is not included in keepedORFs, add it
    if geneID not in keepedORFs.keys():
        keepedORFs[geneID] = {ribORFend: [ribORFstart, nORF, StartCodon] }
    if geneID in keepedORFs.keys():
        # If it's already there check if the same ribORF end is shared
        if ribORFend in keepedORFs[geneID].keys():
            # If it is check if the stating codon is earlier or not
            # if it is store the start and also the ORF number of the gene
            if ATGmode == "F":

                if int(ribORFstart) < int(keepedORFs[geneID][ribORFend][0]):
                    keepedORFs[geneID][ribORFend][0] = ribORFstart
                    keepedORFs[geneID][ribORFend][1] = nORF
                    keepedORFs[geneID][ribORFend][2] = StartCodon
            elif ATGmode == "T":
                if keepedORFs[geneID][ribORFend][2] == "ATG":
                    if int(ribORFstart) < int(keepedORFs[geneID][ribORFend][0]) and StartCodon == "ATG":
                        keepedORFs[geneID][ribORFend][0] = ribORFstart
                        keepedORFs[geneID][ribORFend][1] = nORF
                        keepedORFs[geneID][ribORFend][2] = StartCodon
                else:
                    if int(ribORFstart) < int(keepedORFs[geneID][ribORFend][0]) or StartCodon == "ATG":
                        keepedORFs[geneID][ribORFend][0] = ribORFstart
                        keepedORFs[geneID][ribORFend][1] = nORF
                        keepedORFs[geneID][ribORFend][2] = StartCodon
        else:
            keepedORFs[geneID].update({ribORFend: [ribORFstart, nORF, StartCodon] })
file.close()
file = open(filePath,"r")
for line in file:
    ribORF, chr, strand, transcriptStart, TranscriptEnd, ORFstart, ORFend, exonNumber, exonStarts, exonEnds = line.split("\t")
    geneID, ribORFchr, strand_nORF_transcriptLength, ribORFstart, ribORFend_ORFType_StartCodon = ribORF.split(":")
    ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
    ribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")
    if ORFType not in featureChosenList and "All" not in featureChosenList:
        continue
    if nORF == keepedORFs[geneID][ribORFend][1]:
        outFile.writelines(line)

outFile.close()
file.close()
