# We have to do 2 loops
# THe first one is to annotate all the canonical or prone to be canonical ORFs
# Then if we were not able to find the canonical or prone to be canonical ORF just do not print it.


import sys

try:
    functionName, filePath, outFilePath = sys.argv
except:
    print("error: FixRibORFGenePred.py <Input file> <Output file>")
    quit()

file = open(filePath,"r")
outFile = open(outFilePath,"w+")
canonicalFrame = {}
ORFstartEnd = {}
# function to obtain the right numbers

def realORFcoordinates(exonNumber, strand, transcriptStart, TranscriptEnd, ribORFstart, exonStarts, exonEnds):

    # First of all check for the right strand
    if strand == "+":
        if int(exonNumber) == 1:
            ORFstartNew = int(transcriptStart) + int(ribORFstart) - 1
            ORFendNew = int(transcriptStart) + int(ribORFend) - 1
            return(ORFstartNew, ORFendNew)
        # If has more than one exon
        else:
            exonStartsSplitted = exonStarts.split(",")
            exonEndsSplitted = exonEnds.split(",")
            previousSequences = 0
            ORFstartNew = 0
            ORFendNew = 0
            for i in range(0, (len(exonStartsSplitted) - 1)):
                sizeCurrentBlock = previousSequences + int(exonEndsSplitted[i]) - int(exonStartsSplitted[i])
                if int(ribORFstart) < sizeCurrentBlock and ORFstartNew == 0:
                    ORFstartNew = (int(ribORFstart) - previousSequences) + int(exonStartsSplitted[i]) - 1

                if int(ribORFend) <= (sizeCurrentBlock) and ORFendNew == 0:
                    ORFendNew = (int(ribORFend) - previousSequences) + int(exonStartsSplitted[i]) - 1
                previousSequences = sizeCurrentBlock

            return(ORFstartNew, ORFendNew)
    else:
        if int(exonNumber) == 1:
            ORFendNew = int(TranscriptEnd) - int(ribORFstart) + 1
            ORFstartNew = int(TranscriptEnd) - int(ribORFend) - 2 + 3
            return(ORFstartNew, ORFendNew)
        # Continue from here

        else:
            exonStartsSplitted = exonStarts.split(",")
            exonEndsSplitted = exonEnds.split(",")
            previousSequences = 0
            ORFstartNew = 0
            ORFendNew = 0

            for i in range((len(exonStartsSplitted) - 2), -1, -1):

                sizeCurrentBlock = previousSequences + int(exonEndsSplitted[i]) - int(exonStartsSplitted[i])

                if int(ribORFstart) < sizeCurrentBlock and ORFendNew == 0:

                    ORFendNew = int(exonEndsSplitted[i]) - (int(ribORFstart) - previousSequences) + 1

                if int(ribORFend) <= (sizeCurrentBlock + 1) and ORFstartNew == 0:
                    ORFstartNew = int(exonEndsSplitted[i]) - (int(ribORFend) - previousSequences) + 1
                previousSequences = sizeCurrentBlock
            return (ORFstartNew, ORFendNew)

for line in file:
    ribORF, chr, strand, transcriptStart, TranscriptEnd, ORFstart, ORFend, exonNumber, exonStarts, exonEnds = line.split("\t")
    geneID, ribORFchr, strand_nORF_transcriptLength, ribORFstart, ribORFend_ORFType_StartCodon = ribORF.split(":")
    ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
    ribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")

    ORFtranscriptRelation = (int(ribORFend) - int(ribORFstart)) / int(transcriptLength) * 100

    # Here we fix those transcripts without UTR regions, but only the canonical isoform

    if int(ribORFstart) == 1 and int(transcriptLength) == (int(ribORFend) - 1) and (ORFType == "uORF" or ORFType == "polycistronic"):
        ribORFNEW = geneID +":"+ ribORFchr +":"+ strand_nORF_transcriptLength +":"+ ribORFstart +":"+ ribORFend +"|"+ "canonical" +"|"+ StartCodon
        lineNew = ribORFNEW +"\t"+ chr +"\t"+ strand +"\t"+ transcriptStart +"\t"+ TranscriptEnd +"\t"+ ORFstart +"\t"+ TranscriptEnd +"\t"+ exonNumber +"\t"+ exonStarts +"\t"+ exonEnds
        ORFstartEnd[geneID] = [int(ribORFstart), int(ribORFend)]
        canonicalFrame[geneID] = int(ribORFstart) % 3
        outFile.writelines(lineNew)

    # ORFs without 5'UTR but with 3'UTR annotated
    elif int(ribORFstart) == 1 and ORFtranscriptRelation > 20 and ORFType == "uORF" and strand == "-" and geneID not in canonicalFrame.keys():

        ribORFNEW = geneID +":"+ ribORFchr +":"+ strand_nORF_transcriptLength +":"+ ribORFstart +":"+ ribORFend +"|"+ "canonical" +"|"+ StartCodon
        ORFstartNew, ORFendNew = realORFcoordinates(exonNumber, strand, transcriptStart, TranscriptEnd, ribORFstart, exonStarts, exonEnds)
        lineNew = ribORFNEW +"\t"+ chr +"\t"+ strand +"\t"+ transcriptStart +"\t"+ TranscriptEnd +"\t"+ str(ORFstartNew) + "\t" + str(ORFendNew) +"\t"+ exonNumber +"\t"+ exonStarts +"\t"+ exonEnds
        canonicalFrame[geneID] = int(ribORFstart) % 3
        ORFstartEnd[geneID] = [int(ribORFstart), int(ribORFend)]
        outFile.writelines(lineNew)

    # ORFs without 3'UTR but with 5'UTR annotated
    elif int(transcriptLength) == (int(ribORFend) - 1) and ORFtranscriptRelation > 20 and ORFType == "polycistronic" and strand == "+" and geneID not in canonicalFrame.keys():

        ribORFNEW = geneID + ":" + ribORFchr + ":" + strand_nORF_transcriptLength + ":" + ribORFstart + ":" + ribORFend + "|" + "canonical" + "|" + StartCodon
        ORFstartNew, ORFendNew = realORFcoordinates(exonNumber, strand, transcriptStart, TranscriptEnd, ribORFstart, exonStarts, exonEnds)
        lineNew = ribORFNEW +"\t"+ chr +"\t"+ strand +"\t"+ transcriptStart +"\t"+ TranscriptEnd +"\t"+ str(ORFstartNew) + "\t" + str(ORFendNew) +"\t"+ exonNumber +"\t"+ exonStarts +"\t"+ exonEnds
        canonicalFrame[geneID] = int(ribORFstart) % 3
        ORFstartEnd[geneID] = [int(ribORFstart), int(ribORFend)]
        outFile.writelines(lineNew)

    # Here we fix those transcripts without UTR regions in negative strand that are called by default uORFs

    elif ORFType == "uORF" and strand == "-" and geneID in canonicalFrame.keys():

        if canonicalFrame[geneID] == (int(ribORFstart) % 3) and (int(ribORFstart) >= ORFstartEnd[geneID][0] and int(ribORFend) <= ORFstartEnd[geneID][1]):
            ORFTypeNew = "truncation"
        elif canonicalFrame[geneID] != (int(ribORFstart) % 3) and (int(ribORFstart) >= ORFstartEnd[geneID][0] and int(ribORFend) <= ORFstartEnd[geneID][1]):
            ORFTypeNew = "internal"
        elif int(ribORFstart) < ORFstartEnd[geneID][0] and int(ribORFend) <= ORFstartEnd[geneID][0]:
            ORFTypeNew = "uORF"
        elif int(ribORFstart) < ORFstartEnd[geneID][0] and int(ribORFend) <= ORFstartEnd[geneID][1]:
            ORFTypeNew = "overlap.uORF"
        elif int(ribORFstart) > ORFstartEnd[geneID][1]:
            ORFTypeNew = "polycistronic"
        elif int(ribORFstart) > ORFstartEnd[geneID][0] and int(ribORFend) > ORFstartEnd[geneID][1]:
            ORFTypeNew = "external"
        else:
            ORFTypeNew = "other"
        ribORFNEW = geneID + ":" + ribORFchr + ":" + strand_nORF_transcriptLength + ":" + ribORFstart + ":" + ribORFend + "|" + ORFTypeNew + "|" + StartCodon

        ORFstartNew, ORFendNew = realORFcoordinates(exonNumber, strand, transcriptStart, TranscriptEnd, ribORFstart, exonStarts, exonEnds)
        lineNew = ribORFNEW + "\t" + chr + "\t" + strand + "\t" + transcriptStart + "\t" + TranscriptEnd + "\t" + str(ORFstartNew) + "\t" + str(ORFendNew) + "\t" + exonNumber + "\t" + exonStarts + "\t" + exonEnds
        outFile.writelines(lineNew)

    elif ORFType == "polycistronic" and strand == "+":
        if geneID not in canonicalFrame.keys():
            outFile.writelines(line)
            continue

        if canonicalFrame[geneID] == (int(ribORFstart) % 3):
            ORFTypeNew = "truncation"
        else:
            ORFTypeNew = "internal"

        ribORFNEW = geneID + ":" + ribORFchr + ":" + strand_nORF_transcriptLength + ":" + ribORFstart + ":" + ribORFend + "|" + ORFTypeNew + "|" + StartCodon

        ORFstartNew, ORFendNew = realORFcoordinates(exonNumber, strand, transcriptStart, TranscriptEnd, ribORFstart, exonStarts, exonEnds)
        lineNew = ribORFNEW + "\t" + chr + "\t" + strand + "\t" + transcriptStart + "\t" + TranscriptEnd + "\t" + str(ORFstartNew) + "\t" + str(ORFendNew) + "\t" + exonNumber + "\t" + exonStarts + "\t" + exonEnds
        outFile.writelines(lineNew)
        # Working on this part
    else:
        outFile.writelines(line)
outFile.close()