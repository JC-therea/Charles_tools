import sys

try:
    functionName, refPAth, originalFunPath, filteredFunPath, rescuedOutPath = sys.argv
except:
    print("error: gffCompareExtraction.py <reference.gff> <funannotate.original.gff> <funannotate.filtered.gff> <Output gff>")
    quit()

oldFunGenes = []
oldFunTranscripts = []

# First get the gene and transcript names of the filtered file

gxfOutput = open(rescuedOutPath,"w+")

with open(filteredFunPath, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        if feature == "gene":
            gxfOutput.writelines(line)
            geneName = attr.split("ID=")[1].split(";")[0]
            if geneName not in oldFunGenes:
                oldFunGenes.append(geneName)
        elif feature in ["mRNA", "rRNA", "ncRNA", "tRNA"]:
            gxfOutput.writelines(line)
            transcriptName = attr.split("ID=")[1].split(";")[0]
            geneName = attr.split("Parent=")[1].split(";")[0]
            if transcriptName not in oldFunTranscripts:
                oldFunTranscripts.append(transcriptName)

        else:
            gxfOutput.writelines(line)
            transcriptName = attr.split("Parent=")[1].split(";")[0]
            if transcriptName not in oldFunTranscripts:
                oldFunTranscripts.append(transcriptName)

utr3Dict = {}
utr5Dict = {}

Transcript2Gene = {}

# From the original Funannotate file (without the filtering)
# create a dictionary to recover the gene name from the transcript name

with open(originalFunPath, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        if feature == "mRNA":
            transcriptName = attr.split("ID=")[1].split(";")[0]
            geneName = attr.split("Parent=")[1].split(";")[0]
            Transcript2Gene[transcriptName] = geneName

# Store UTRs of all the transcripts
# From the original file
with open(originalFunPath, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        if feature == "five_prime_UTR":
            transcript = attr.split("Parent=")[1].split(";")[0]
            gene = Transcript2Gene[transcript]
            # Keep the longest 5' UTR
            if gene in utr5Dict.keys():
                if int(start) < int(utr5Dict[gene][0]):
                    utr5Dict[gene][0] = start
                if int(end) > int(utr5Dict[gene][1]):
                    utr5Dict[gene][1] = end
            else:
                utr5Dict[gene] = [start, end]

        elif feature == "three_prime_UTR":
            transcript = attr.split("Parent=")[1].split(";")[0]
            gene = Transcript2Gene[transcript]
            # Keep the longest 3' UTR
            if gene in utr3Dict.keys():
                if int(start) < int(utr3Dict[gene][0]):
                    utr3Dict[gene][0] = start
                if int(end) > int(utr3Dict[gene][1]):
                    utr3Dict[gene][1] = end
            else:
                utr3Dict[gene] = [start, end]

# Keep all the genes and transcript of the reference sequence

exonStorage = {}
transcripts_with_CDS = []
genes_with_CDS = []
Transcript2geneRef = {}

with open(refPAth, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")

        if feature == "mRNA":
            refGene = attr.split("Parent=")[1]
            if ";" in refGene:
                refGene = refGene.split(";")[0]
            else:
                refGene = refGene.split("\n")[0]
            refTranscript = attr.split("ID=")[1].split(";")[0]

            if refTranscript not in Transcript2geneRef.keys():
                Transcript2geneRef[refTranscript] = refGene

# Store exons because we only want to modify the ends of the transcript
# and if they are multiexonic we only want to modify one exon
# Also keep the transcripts and genes that have CDS because we want those

with open(refPAth, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")

        if feature in ["exon"]:
            refTranscript = attr.split("Parent=")[1]
            if ";" in refTranscript:
                refTranscript = refTranscript.split(";")[0]
            else:
                refTranscript = refTranscript.split("\n")[0]
            if refTranscript in Transcript2geneRef.keys():
                refGene = Transcript2geneRef[refTranscript]
                if refGene not in oldFunGenes:
                    if refGene not in exonStorage.keys():
                        exonStorage[refGene] = [start,end]
                    else:
                        exonStorage[refGene].append(start)
                        exonStorage[refGene].append(end)

        elif feature in ["CDS"]:
            refTranscript = attr.split("Parent=")[1]
            if ";" in refTranscript:
                refTranscript = refTranscript.split(";")[0]
            else:
                refTranscript = refTranscript.split("\n")[0]
            if refTranscript not in transcripts_with_CDS:
                transcripts_with_CDS.append(refTranscript)
                if refTranscript in Transcript2geneRef.keys():
                    genes_with_CDS.append(Transcript2geneRef[refTranscript])

exonCounting = {}

# Here is the main part where we rescue the transcripts

with open(refPAth, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")

        # Rescue the gene and modify it if they have UTRs

        if feature == "gene":
            refGene = attr.split("ID=")[1]
            if ";" in refGene:
                refGene = refGene.split(";")[0]
            else:
                refGene = refGene.split("\n")[0]
            if refGene not in oldFunGenes and refGene in genes_with_CDS:
                if strand == "+":
                    if refGene in utr5Dict.keys():
                        if int(start) > int(utr5Dict[refGene][0]):
                            start = utr5Dict[refGene][0]
                    if refGene in utr3Dict.keys():

                        if int(end) < int(utr3Dict[refGene][1]):
                            end = utr3Dict[refGene][1]
                else:
                    if refGene in utr3Dict.keys():
                        if int(start) > int(utr3Dict[refGene][0]):
                            start = utr3Dict[refGene][0]
                    if refGene in utr5Dict.keys():
                        if int(end) < int(utr5Dict[refGene][1]):
                            end = utr5Dict[refGene][1]
                new_attr = "ID=" + refGene + ";" + "\n"
                new_line = "\t".join(map(str,[chr, source, feature, start, end, score, strand, frame, new_attr]))
                gxfOutput.writelines(new_line)

        # Rescue mRNAs and modify them if they have UTRs

        elif feature == "mRNA":
            refGene = attr.split("Parent=")[1]
            if ";" in refGene:
                refGene = refGene.split(";")[0]
            else:
                refGene = refGene.split("\n")[0]
            refTranscript = attr.split("ID=")[1].split(";")[0]
            if refGene not in oldFunGenes and refTranscript in transcripts_with_CDS:
                if strand == "+":
                    if refGene in utr5Dict.keys():
                        if int(start) > int(utr5Dict[refGene][0]):
                            start = utr5Dict[refGene][0]
                    if refGene in utr3Dict.keys():
                        if int(end) < int(utr3Dict[refGene][1]):
                            end = utr3Dict[refGene][1]
                else:
                    if refGene in utr3Dict.keys():
                        if int(start) > int(utr3Dict[refGene][0]):
                            start = utr3Dict[refGene][0]
                    if refGene in utr5Dict.keys():
                        if int(end) < int(utr5Dict[refGene][1]):
                            end = utr5Dict[refGene][1]
                new_attr = "ID=" + refGene + "-T1;Parent=" + refGene + ";product=hypothetical protein;" + "\n"
                new_line = "\t".join(map(str,[chr, source, feature, start, end, score, strand, frame, new_attr]))
                gxfOutput.writelines(new_line)

        # Rescue the exons and modify them if they have UTRs

        elif feature == "exon":
            refTranscript = attr.split("Parent=")[1]
            if ";" in refTranscript:
                refTranscript = refTranscript.split(";")[0]
            else:
                refTranscript = refTranscript.split("\n")[0]
            if refTranscript in Transcript2geneRef.keys():
                refGene = Transcript2geneRef[refTranscript]
            else:
                continue
            if (refGene not in oldFunGenes) and (refTranscript in transcripts_with_CDS) and (refTranscript in Transcript2geneRef.keys()):
                # refGene = Transcript2geneRef[refTranscript]
                if len(exonStorage[refGene]) == 2:
                    if strand == "+":
                        if refGene in utr5Dict.keys():
                            if int(start) > int(utr5Dict[refGene][0]):
                                start = utr5Dict[refGene][0]
                        if refGene in utr3Dict.keys():
                            if int(end) < int(utr3Dict[refGene][1]):
                                end = utr3Dict[refGene][1]
                    else:
                        if refGene in utr3Dict.keys():
                            if int(start) > int(utr3Dict[refGene][0]):
                                start = utr3Dict[refGene][0]
                        if refGene in utr5Dict.keys():
                            if int(end) < int(utr5Dict[refGene][1]):
                                end = utr5Dict[refGene][1]
                else:
                    if strand == "+":
                        if start == min(exonStorage[refGene]):
                            if refGene in utr5Dict.keys():
                                if int(start) > int(utr5Dict[refGene][0]):
                                    start = utr5Dict[refGene][0]
                        if end == max(exonStorage[refGene]):
                            if refGene in utr3Dict.keys():
                                if int(end) < int(utr3Dict[refGene][1]):
                                    end = utr3Dict[refGene][1]
                    else:
                        if start == min(exonStorage[refGene]):
                            if refGene in utr3Dict.keys():
                                if int(start) > int(utr3Dict[refGene][0]):
                                    start = utr3Dict[refGene][0]
                        if end == max(exonStorage[refGene]):
                            if refGene in utr5Dict.keys():
                                if int(end) < int(utr5Dict[refGene][1]):
                                    end = utr5Dict[refGene][1]
                if refGene not in exonCounting.keys():
                    exonCounting[refGene] = 1
                else:
                    exonCounting[refGene] += 1
                new_attr = "ID=" + refGene + "-T1.exon" + str(exonCounting[refGene]) + ";Parent=" + refGene + "-T1;" + "\n"
                new_line = "\t".join(map(str, [chr, source, feature, start, end, score, strand, frame, new_attr]))
                gxfOutput.writelines(new_line)

        # Get the CDS and keep it as it is

        elif feature in ["CDS"]:
            refTranscript = attr.split("Parent=")[1]
            if ";" in refTranscript:
                refTranscript = refTranscript.split(";")[0]
            else:
                refTranscript = refTranscript.split("\n")[0]
            if refGene not in oldFunGenes and refTranscript in Transcript2geneRef.keys():
                refGene = Transcript2geneRef[refTranscript]
                new_attr = "ID=" + refGene + "-T1.cds;Parent=" + refGene + "-T1;" + "\n"
                new_line = "\t".join(map(str, [chr, source, feature, start, end, score, strand, frame, new_attr]))
                gxfOutput.writelines(new_line)


# Now print UTRs with agat

gxfOutput.close()
# Check which transcripts are stored in the file
