
# With the help of ChatGPT
# read the input file into memory

import sys

try:
    functionName, refPAth, rescuedOutPath, Chopped_UTRs = sys.argv
except:
    print("error: ChopityChop.py <reference.gff> <funannotate.rescued.gff> <Output gff>")
    quit()


Ori = refPAth
FunRescue = rescuedOutPath
with open(Ori,  "r") as input_file:
    lines = input_file.readlines()

# initialize a dictionary to store the CDS and tRNA features
cds_trna_dict = {}
transcripts2genes = {}
genes2transcript = {}
transcripts_with_UTRs = []
# loop over each line in the input file and store the CDS and tRNA features
for line in lines:
    line = line.strip()

    # skip comment lines
    if line.startswith("#"):
        continue

    # split the line into fields
    chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")

    # check if the feature is a CDS or tRNA
    if feature in ["tRNA"]:
        # add the start and end positions to the dictionary
        if (chr, strand) in cds_trna_dict:
            cds_trna_dict[(chr, strand)].append((int(start), int(end)))
        else:
            cds_trna_dict[(chr, strand)] = [(int(start), int(end))]

with open(FunRescue,  "r") as input_file:
    lines = input_file.readlines()

    for line in lines:
        line = line.strip()

        # skip comment lines
        if line.startswith("#"):
            continue

        # split the line into fields
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        fields = line.split("\t")

        if feature in ["three_prime_UTR","five_prime_UTR"]:
            transcript = attr.split("Parent=")[1].split(";")[0]
            transcripts_with_UTRs.append(transcript)
        if feature in ["tRNA", "mRNA"]:
            transcriptName = attr.split("ID=")[1].split(";")[0]
            geneName = attr.split("Parent=")[1].split(";")[0]
            transcripts2genes[transcriptName] = geneName
            genes2transcript[geneName] = transcriptName

# open the output file
with open(Chopped_UTRs, "w") as output_file:
    # loop over each line in the input file
    for line in lines:
        line = line.strip()

        # skip comment lines
        if line.startswith("#"):
            output_file.write(line + "\n")
            continue

        # split the line into fields
        chr, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        fields = line.split("\t")
        # check if the feature is a UTR, gene, or mRNA
        if feature in ["three_prime_UTR","five_prime_UTR", "gene", "mRNA", "exon"]:

            if feature == "gene":
                geneName = attr.split("ID=")[1].split(";")[0]
                if geneName in genes2transcript.keys():
                    transcriptName = genes2transcript[geneName]
                else:
                    continue

            elif feature == "mRNA":
                transcriptName = attr.split("ID=")[1].split(";")[0]
            elif feature in ["three_prime_UTR","five_prime_UTR", "exon"]:
                transcriptName = attr.split("Parent=")[1].split(";")[0]
            # loop over the CDS and tRNA features on the same chromosome and strand
            if transcriptName not in transcripts_with_UTRs:
                output_file.write(line + "\n")
                continue
            for (cds_start, cds_end) in cds_trna_dict.get((chr, strand), []):
                # check if the features overlap
                if int(start) == cds_start and int(end) == cds_end:
                    continue
                # Avoid modifying annotated tRNAs
                elif attr.startswith("ID=gene-YNC"):
                    continue
                elif int(start) < cds_end and int(end) > cds_start:
                    print(attr)
                    # update the start and end positions of the feature

                    if strand == "+":
                        if int(start) < cds_start and int(end) > cds_end and strand == "+":
                            leftRest = cds_start - int(start)
                            rightRest = int(end) - cds_end
                            if feature == "five_prime_UTR":
                                start = cds_end + 1
                            elif feature == "three_prime_UTR":
                                end = cds_start - 1
                            elif leftRest < rightRest:
                                start = cds_end + 1
                            else:
                                end = cds_start - 1
                        else:
                            start = max(int(start), cds_end) + 1
                    elif strand == "-":
                        if int(start) < cds_start and int(end) > cds_end and strand == "+":
                            leftRest = cds_start - int(start)
                            rightRest = int(end) - cds_end
                            if feature == "five_prime_UTR":
                                start = cds_start - 1
                            elif feature == "three_prime_UTR":
                                end = cds_end + 1
                            elif leftRest < rightRest:
                                start = cds_end + 1
                            else:
                                end = cds_start - 1
                        else:
                            end = min(int(end), cds_start) - 1
            # write out the modified feature
            output_file.write("\t".join(fields[:3] + [str(start), str(end)] + fields[5:]) + "\n")

        # if the feature is not a UTR, gene, or mRNA, write it out as is
        else:
            output_file.write(line + "\n")
