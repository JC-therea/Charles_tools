import sys

try:
    functionName, gffCompareTracking_path, symbolsChosen, gxfToExtract_path, gxfOutput_path = sys.argv
except:
    print("error: gffCompareExtraction.py <gffCompare.tracking file> <gffCompare identifiers chosen separated by ','> <gtf/gff to extract sequences> <Output gff/gtf>")
    quit()

if "," in symbolsChosen:
    symbolsChosen_list = symbolsChosen.split(",")
else:
    symbolsChosen_list = [symbolsChosen]

desiredTranscripts = []

with open(gffCompareTracking_path, "r") as file:
    for line in file:
        if line.split("\t")[3] in symbolsChosen_list:
            desiredTranscripts.append(line.split("\t")[4].split(":")[1].split("|")[0])
            desiredTranscripts.append(line.split("\t")[4].split(":")[1].split("|")[1])

# Depend on extension

gxfOutput = open(gxfOutput_path, "w+")

with open(gxfToExtract_path, "r") as file:
    for line in file:
        if "#" in line[0]:
            continue
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        for transcript in desiredTranscripts:
            if transcript in attr:
                gxfOutput.writelines(line)
                break
        #if gxfToExtract_path.split(".")[-1] == "gtf":

        #elif gxfToExtract_path.split(".")[-1] == "gff" or gxfToExtract_path.split(".")[-1] == "gff3":
gxfOutput.close()
