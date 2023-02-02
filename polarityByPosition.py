from Bio import SeqIO
import pandas as pd
import math
import argparse

parser = argparse.ArgumentParser(description="This program calculates the polarity distribution of the amino acids",
								 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input", default="", type=str, help="path of the input fasta file", required=True)
parser.add_argument("-w", "--windows", default=10, type=int, help="Set the number of windows that you analyze, if you have smaller peptides than the windows probably there is going to be NaNs", required=True)
parser.add_argument("-o", "--output", default="", type=str, help="Name of the output file", required=False)
args = parser.parse_args()
inputPath = args.input
outPath = args.output
NPoints = args.windows

# inputPath = "/home/jmontanes/Documents/IQtree_Gene_duplication/Insects/Standard_annotation/OutputsR/AA_content/AA_alignments/PolarityDist/Dmel_N0.fa"
# outPath = ""

fasta_sequences = SeqIO.parse(open(inputPath),'fasta')
records = len([1 for line in open(inputPath) if line.startswith(">")])
aa_type={
"A" : "nonPolar", "C" : "Polar",  "D" : "Acidic", "E" : "Acidic", "F" : "nonPolar", "G" : "nonPolar", "H" : "Basic", "I" : "nonPolar", "K" : "Basic", "L" : "nonPolar",
"M" : "nonPolar", "N" : "Polar", "P" : "nonPolar", "Q" : "Polar", "R" : "Basic", "S" : "Polar", "T" : "Polar", "V" : "nonPolar", "W" : "nonPolar", "Y" : "Polar"
}

PolarList = ["C","N","Q","S","T","Y"]
nonPolarList = ["A","F","G","I","L","M","P","V","W"]
AcidicList = ["D","E"]
BasicList = ["H","K","R"]


dfPolar = pd.DataFrame(index=range(records),columns=range(NPoints))
dfNonPolar = pd.DataFrame(index=range(records),columns=range(NPoints))
dfAcidic = pd.DataFrame(index=range(records),columns=range(NPoints))
dfBasic = pd.DataFrame(index=range(records),columns=range(NPoints))

SeqNumber = 0
fastaNames = []
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    #print(name)
    if "." in sequence:
        print("Not allowed symbols")
        break
    fastaNames.append(name)
    RelStart = 0
    WholeRowPolar = {}
    WholeRowNonPolar = {}
    WholeRowAcidic = {}
    WholeRowBasic = {}

    for relPos in range(NPoints):
        #print(relPos)
        AbsPosStart = RelStart
        if relPos != (NPoints - 1):
            AbsPosEnd = math.floor((relPos + 1) / NPoints * len(sequence))
        else:
            AbsPosEnd = len(sequence) - 1

        Subseq = sequence[AbsPosStart:AbsPosEnd]
        Polar = 0
        Acidic = 0
        Basic = 0
        NonPolar = 0

        for letter in Subseq:
            if letter in PolarList:
                Polar += 1
            elif letter in AcidicList:
                Acidic += 1
            elif letter in BasicList:
                Basic += 1
            elif letter in nonPolarList:
                NonPolar += 1

        totAAs = Polar + NonPolar + Acidic + Basic
        if totAAs == 0:
            continue

        WholeRowPolar[relPos] = Polar/totAAs
        WholeRowNonPolar[relPos] = NonPolar / totAAs
        WholeRowAcidic[relPos] = Acidic / totAAs
        WholeRowBasic[relPos] = Basic / totAAs

        RelStart = AbsPosEnd

    dfPolar.loc[SeqNumber] = pd.Series(WholeRowPolar)
    dfAcidic.loc[SeqNumber] = pd.Series(WholeRowAcidic)
    dfBasic.loc[SeqNumber] = pd.Series(WholeRowBasic)
    dfNonPolar.loc[SeqNumber] = pd.Series(WholeRowNonPolar)
    SeqNumber += 1


dfPolarDef = dfPolar.set_index(pd.Series(fastaNames))
dfPolarDef.to_csv(path_or_buf= outPath + "Polar.tsv", sep = "\t")

dfAcidicDef = dfAcidic.set_index(pd.Series(fastaNames))
dfAcidicDef.to_csv(path_or_buf= outPath + "Acidic.tsv", sep = "\t")

dfBasicDef = dfBasic.set_index(pd.Series(fastaNames))
dfBasicDef.to_csv(path_or_buf= outPath + "Basic.tsv", sep = "\t")

dfNonPolarDef = dfNonPolar.set_index(pd.Series(fastaNames))
dfNonPolarDef.to_csv(path_or_buf= outPath + "nonPolar.tsv", sep = "\t")


