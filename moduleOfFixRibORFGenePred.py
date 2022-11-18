import sys

try:
    functionName, filePath, outFilePath = sys.argv
except:
    print("error: moduleOfFixRibORFGenePred.py <Input file> <Output file>")
    quit()

file = open(filePath,"r")
outFile = open(outFilePath,"w+")
canonicalFrame = {}
ORFstartEnd = {}

keepedORFs = {}
#allowed_ORFtypes = ["canonical", "external", "internal", "overlap.uORF", "polycistronic", "uORF", "noncoding", "extension"]

for line in file:
    ribORF, chr, strand, transcriptStart, TranscriptEnd, ORFstart, ORFend, exonNumber, exonStarts, exonEnds = line.split("\t")
    geneID, ribORFchr, strand_nORF_transcriptLength, ribORFstart, ribORFend_ORFType_StartCodon = ribORF.split(":")
    ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
    ribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")
    #if ORFType not in allowed_ORFtypes:
    #    continue
    if geneID not in keepedORFs.keys():
        keepedORFs[geneID] = {ribORFend: [ribORFstart,nORF] }
    if geneID in keepedORFs.keys():
        if ribORFend in keepedORFs[geneID].keys():
            if int(ribORFstart) < int(keepedORFs[geneID][ribORFend][0]):
                keepedORFs[geneID][ribORFend][0] = ribORFstart
                keepedORFs[geneID][ribORFend][1] = nORF
        else:
            keepedORFs[geneID].update({ribORFend: [ribORFstart,nORF] })
file.close()
file = open(filePath,"r")
for line in file:
    ribORF, chr, strand, transcriptStart, TranscriptEnd, ORFstart, ORFend, exonNumber, exonStarts, exonEnds = line.split("\t")
    geneID, ribORFchr, strand_nORF_transcriptLength, ribORFstart, ribORFend_ORFType_StartCodon = ribORF.split(":")
    ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
    ribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")
    #print(ribORF)
    if nORF == keepedORFs[geneID][ribORFend][1]:
        outFile.writelines(line)

file.close()
