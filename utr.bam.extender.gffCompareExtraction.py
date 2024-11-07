import sys

try:
    functionName, gffCompareTracking_path, symbolsChosen, gxfToExtract_path, gxfOutput_path = sys.argv
except:
    print("error: gffCompareExtraction.py <gffCompare.tracking file> <gffCompare identifiers chosen separated by ','> <gff to extract sequences> <Output gff>")
    quit()

if "," in symbolsChosen:
    symbolsChosen_list = symbolsChosen.split(",")
else:
    symbolsChosen_list = [symbolsChosen]

desiredTranscripts = []
desiredGenes = []
with open(gffCompareTracking_path, "r") as file:
    for line in file:
        if line.split("\t")[3] in symbolsChosen_list:
            desiredGenes.append(line.split("\t")[4].split(":")[1].split("|")[0])
            desiredTranscripts.append(line.split("\t")[4].split(":")[1].split("|")[1])

geneTranscripts = {}
multiTranscriptsValids = []
validTranscripts = []
exonsInTranscript = {}

# Check if there is more than one transcript per accepted gene
with open(gxfToExtract_path, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        if feature in ["mRNA", "rRNA", "ncRNA", "tRNA","telomerase_RNA"]:
            transcriptName = attr.split("ID=")[1].split(";")[0]
            geneName = attr.split("Parent=")[1].split(";")[0]
            if transcriptName in desiredTranscripts:
                if geneName not in geneTranscripts.keys():
                    geneTranscripts[geneName] = [transcriptName]
                else:
                    geneTranscripts[geneName].append(transcriptName)
                    multiTranscriptsValids.append(geneName)

# Store the exons per transcripts only of multitranscript genes
with open(gxfToExtract_path, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        if feature == "exon":
            transcriptName = attr.split("Parent=")[1].split(";")[0]
            if transcriptName in desiredTranscripts:
                if transcriptName not in exonsInTranscript.keys():
                    exonsInTranscript[transcriptName] = 1
                else:
                    exonsInTranscript[transcriptName] += 1

# Pick one transcript (the one with less exons) per gene in problematic genes
for gene in geneTranscripts.keys():
    if gene in multiTranscriptsValids:
        numExons = 999
        for transcript in geneTranscripts[gene]:
            if exonsInTranscript[transcript] < numExons:
                goodTranscript = transcript
                numExons = exonsInTranscript[transcript]
        validTranscripts.append(goodTranscript)

    else:
        validTranscripts.append(geneTranscripts[gene][0])

gxfOutput = open(gxfOutput_path, "w+")

with open(gxfToExtract_path, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        if feature == "gene":
            geneName = attr.split("ID=")[1].split(";")[0]
            if geneName in desiredGenes:
                gxfOutput.writelines(line)
        elif feature in ["mRNA", "rRNA", "ncRNA", "tRNA","telomerase_RNA"]:
            transcriptName = attr.split("ID=")[1].split(";")[0]
            if transcriptName in desiredTranscripts and transcriptName in validTranscripts:
                gxfOutput.writelines(line)
        else:
            transcriptName = attr.split("Parent=")[1].split(";")[0]
            if transcriptName in desiredTranscripts and transcriptName in validTranscripts:
                gxfOutput.writelines(line)
gxfOutput.close()
