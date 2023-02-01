from Bio import SeqIO
import pandas as pd
import math
import argparse

parser = argparse.ArgumentParser(description="This program calculates the polarity distribution of the amino acids",
								 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input", default="", type=str, help="path of the input fasta file", required=True)
parser.add_argument("-o", "--output", default="", type=str, help="Name of the output file", required=False)

args = parser.parse_args()
inputPath = args.input
fasta_sequences = SeqIO.parse(open(inputPath),'fasta')
records = len([1 for line in open(inputPath) if line.startswith(">")])
outPath = args.output
aa_type={
"A" : "nonPolar", "C" : "Polar",  "D" : "Acidic", "E" : "Acidic", "F" : "nonPolar", "G" : "nonPolar", "H" : "Basic", "I" : "nonPolar", "K" : "Basic", "L" : "nonPolar",
"M" : "nonPolar", "N" : "Polar", "P" : "nonPolar", "Q" : "Polar", "R" : "Basic", "S" : "Polar", "T" : "Polar", "V" : "nonPolar", "W" : "nonPolar", "Y" : "Polar"
}

PolarList = ["C","N","Q","S","T","Y"]
nonPolarList = ["A","F","G","I","L","M","P","V","W"]
AcidicList = ["D","E"]
BasicList = ["H","K","R"]

dfPolar = pd.DataFrame(index=range(records),columns=range(100))
dfNonPolar = pd.DataFrame(index=range(records),columns=range(100))
dfAcidic = pd.DataFrame(index=range(records),columns=range(100))
dfBasic = pd.DataFrame(index=range(records),columns=range(100))

SeqNumber = 0
fastaNames = []
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    print(name)
    if "." in sequence:
        print("Not allowed symbols")
        break
    fastaNames.append(name)
    RelStart = 0
    for relPos in range(100):
        print(relPos)
        AbsPosStart = RelStart
        if relPos != 99:
            AbsPosEnd = math.floor((relPos + 1) / 100 * len(sequence))
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
        dfPolar.loc[SeqNumber,RelStart] = Polar/totAAs
        dfAcidic.loc[SeqNumber, RelStart] = Acidic / totAAs
        dfBasic.loc[SeqNumber, RelStart] = Basic / totAAs
        dfNonPolar.loc[SeqNumber, RelStart] = NonPolar / totAAs

        RelStart = AbsPosEnd

    SeqNumber += 1

dfPolar.set_index(fastaNames)
dfPolar.to_csv(path_or_buf= outPath + "Polar.csv", sep = "\t")

dfAcidic.set_index(fastaNames)
dfAcidic.to_csv(path_or_buf= outPath + "Acidic.csv", sep = "\t")

dfBasic.set_index(fastaNames)
dfBasic.to_csv(path_or_buf= outPath + "Basic.csv", sep = "\t")

dfNonPolar.set_index(fastaNames)
dfNonPolar.to_csv(path_or_buf= outPath + "nonPolar.csv", sep = "\t")


