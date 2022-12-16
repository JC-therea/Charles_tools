from Bio import AlignIO
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="This programs needs the config file to work with the desired species",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c", "--configFile", default="0", type=str, help="path of the configuration file")
parser.add_argument("-m", "--mode", default=0, type=int, help="mode to work with, can be used for alignments among different species (0) or duplications of the same species (1) ")
parser.add_argument("-o", "--outFile", default="", type=str, help="Name of the output as a tsv file")

args = parser.parse_args()

try:
    configFileOp=open(args.configFile, "r")
except:
    print("Not config file found, creating file...")
    configFileOp=open("configFile.txt", "w+")
    configFileOp.writelines("Nodes=")
    configFileOp.writelines("Species=")
    configFileOp.writelines("MSAdir=")
    configFileOp.writelines("SummaryTable=")
    configFileOp.close()
    quit()
try:
    configFileOp=open(args.configFile, "r").readlines()
    attributes = []
    for line in configFileOp:
        attributes.append(line.split("\n")[0])
    Nodes = attributes[0].split("=")[1].split(",")
    Species = attributes[1].split("=")[1].split(",")
    MSAdir = attributes[2].split("=")[1]
    SummaryTable = attributes[3].split("=")[1]
except:
    print("Bad configuration of the config file, unable to run")
    quit()

args = parser.parse_args()
#config = vars(args)
#print(config)

# To do so we will extract the information from a file that I created in R

summaryTablePath = SummaryTable
summaryTable = pd.read_csv(summaryTablePath, sep="\t")

#Nodes = ["N1", "Conserved"]
if args.mode == 0:
    print("Mode comparison between 2 species activated")
    sps1 = Species[0]
    OGs_used = 0
    try:
        sps2 = Species[1]
    except:
        print("Only one species provided, please select 2 species for the analysis")
        quit()

    for node in Nodes:
        OGN1 = summaryTable[summaryTable.AgeNode == node].OrthoGroup
        #MSAdir = "Standard_annotation/Proteomes/OrthoFinder/Results_MSA_IQtree_bigmem/MultipleSequenceAlignments/"

        AlignmentsPais = {}

        for OG in OGN1:
            filePath = MSAdir + OG + ".fa"
            try:
                alignment = AlignIO.read(open(filePath), "fasta")
            except:
                continue
            algn_len = alignment.get_alignment_length()

            SpsProts1 = 0
            SpsProts2 = 0

            # Check the number of proteins per species
            # If there is a different number of proteins than one then avoid the orthogroup

            for sp in range(0, len(alignment)):

                if sps1 in alignment[sp].id:
                    SpsProts1 += 1
                if sps2 in alignment[sp].id:
                    SpsProts2 += 1

            if SpsProts1 != 1 or SpsProts2 != 1:
                # print(filePath + " does not satisfies the condition")
                continue
            OGs_used += 1
            # print(OG)
            for i in range(0, algn_len):
                lettersAlgn = ""
                for sp in range(0, len(alignment)):

                    if sps1 in alignment[sp].id or sps2 in alignment[sp].id:
                        lettersAlgn = lettersAlgn + alignment[sp].seq[i]

                # If all the letters are the same it is a perfect alignment
                if len("".join(set(lettersAlgn))) == 1:
                    AlgnComb = 'Perf'

                else:
                    AlgnComb = ''.join(sorted(lettersAlgn))
                # If the alignment includes missmatches do not save them
                if "-" in AlgnComb:
                    continue
                elif AlgnComb in AlignmentsPais.keys():
                    AlignmentsPais[AlgnComb] += 1
                else:
                    AlignmentsPais[AlgnComb] = 1

        outFilepd = pd.DataFrame(AlignmentsPais.items(), columns=['Alignment', 'Hits'])
        outFilePath = args.outFile + node + "_" + sps1 + "_" + sps2 + ".tsv"
        print(f' For  node {node} {OGs_used} orthogroups were used')
        outFilepd.to_csv(outFilePath, sep="\t", index=False)

if args.mode == 1:
    print("Mode comparison between duplicates activated")
    sps1 = Species[0]
    OGs_used = 0
    for node in Nodes:
        OGN1 = summaryTable[summaryTable.DuplNode == node].OrthoGroup

        AlignmentsPais = {}

        for OG in OGN1:
            filePath = MSAdir + OG + ".fa"
            try:
                alignment = AlignIO.read(open(filePath), "fasta")
            except:
                continue
            algn_len = alignment.get_alignment_length()

            SpsProts1 = 0

            # Check the number of proteins per species
            # If there is a different number of proteins than two then avoid the orthogroup

            for sp in range(0, len(alignment)):

                if sps1 in alignment[sp].id:
                    SpsProts1 += 1

            if SpsProts1 != 2:
                # print(filePath + " does not satisfies the condition")
                continue
            OGs_used += 1
            # print(OG)
            for i in range(0, algn_len):
                lettersAlgn = ""
                for sp in range(0, len(alignment)):

                    if sps1 in alignment[sp].id:
                        lettersAlgn = lettersAlgn + alignment[sp].seq[i]

                # If all the letters are the same it is a perfect alignment
                if len("".join(set(lettersAlgn))) == 1:
                    AlgnComb = 'Perf'

                else:
                    AlgnComb = ''.join(sorted(lettersAlgn))
                # If the alignment includes missmatches do not save them
                if "-" in AlgnComb:
                    continue
                elif AlgnComb in AlignmentsPais.keys():
                    AlignmentsPais[AlgnComb] += 1
                else:
                    AlignmentsPais[AlgnComb] = 1

        outFilepd = pd.DataFrame(AlignmentsPais.items(), columns=['Alignment', 'Hits'])
        outFilePath = args.outFile + node + "_" + sps1  + "_singleDuplicates.tsv"
        print(f' For node {node} {OGs_used} orthogroups were used')
        outFilepd.to_csv(outFilePath, sep="\t", index=False)