
# Run on pythonPackages environment

import sys
import pandas as pd
import warnings
warnings.filterwarnings("ignore")


try:
    functionName, MGCAT_file_path, originalgff_path, gffOut_path = sys.argv
except:
    print("error: AnnotationQC.py <MGCAT file> <original gff> <Output file>")
    quit()


chrList = []

with open(originalgff_path, "r") as file:
   for line in file:
        Chr, Source, Feature, Start, End, Score, Strand, Frame, Attributes = line.split("\t")
        if Chr not in chrList:
            chrList.append(Chr)

chrList.sort()
StoredBlocks = {}
blockNum = 0
gffOut = open(gffOut_path, "w+")
with open(MGCAT_file_path, "r") as file:
   for line in file:
        try:
           chr = int(line.split("\t")[0])
        except:
            continue

        MGCATchr = line.split("\t")[9]
        MGCATstart = line.split("\t")[12]
        MGCATend = line.split("\t")[13]
        MGCATchrAdapted = chrList[int(MGCATchr) - 1]

        if MGCATchrAdapted in StoredBlocks.keys():
            if MGCATstart not in StoredBlocks[MGCATchrAdapted]:
                print(MGCATchrAdapted)
                print(MGCATchr)
                blockNum += 1
                StoredBlocks[MGCATchrAdapted].append(MGCATstart)
                attribute_out = "ID=Syntenic_block_num_"+ str(blockNum) + ";locus=Syntenic_block_num_"+ str(blockNum)
                outLine = MGCATchrAdapted + '\t' + "MGCAT" + '\t' + "SyntenicBlock" + '\t' + MGCATstart + '\t' + MGCATend + '\t' + '.' + '\t' + "+" + '\t' + "." + '\t' + attribute_out + "\n"
                gffOut.writelines(outLine)
        else:
            blockNum += 1
            StoredBlocks[MGCATchrAdapted] = [MGCATstart]
            attribute_out = "ID=Syntenic_block_num_" + str(blockNum) + ";locus=Syntenic_block_num_" + str(blockNum)
            outLine = MGCATchrAdapted + '\t' + "MGCAT" + '\t' + "SyntenicBlock" + '\t' + MGCATstart + '\t' + MGCATend + '\t' + '.' + '\t' + "+" + '\t' + "." + '\t' + attribute_out + "\n"
            gffOut.writelines(outLine)
gffOut.close()