import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="This programs needs the config file to work with the desired species",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--orfology", default="", type=str, help="path to the output tsv file from Conservation_ORFology_expandedFromCandidates")
parser.add_argument("-msa", "--MSAdir", default="", type=str, help="path to the directory where are all the regions aligned")
parser.add_argument("-o", "--outFile", default="", type=str, help="Name of the output file")

args = parser.parse_args()

orfology_table_path = args.orfology
MSA_dir_path = args.MSAdir

# orfology_table_path = "/home/jmontanes/Documents/EvolutionNanopore/Outputs/EvolutionNanopore/Outputs/Transcripts/ConservedORFs_full.tsv"
# MSA_dir_path = "/home/jmontanes/Documents/EvolutionNanopore/Outputs/EvolutionNanopore/Outputs/Transcripts/transcriptAlignments/Fasta"

if not MSA_dir_path.endswith("/"):
    MSA_dir_path = MSA_dir_path + "/"

outPath = args.outFile

########################################################################################
################################ Functions #############################################
########################################################################################

### To extract positions in MSA

def GetPerIdPerGap(refSeq,outSeq):
    if len(refSeq) != len(outSeq):
        print("Different length!")
        exit
    perfMatch = 0.0
    missMatch = 0.0
    gapRef = 0.0
    gapOut = 0.0
    lenAlignment = len(refSeq)

    # Check position by position the identity

    for i in range(0, len(refSeq)):
        refLetter = refSeq[i]
        outLetter = outSeq[i]

        if refLetter == outLetter:
            perfMatch += 1
        else:
            if refLetter == "-" and outLetter != "-":
                gapRef += 1
            elif refLetter != "-" and outLetter == "-":
                gapOut += 1
            elif refLetter != "-" and outLetter != "-":
                missMatch += 1
            elif refLetter == "-" and outLetter == "-":
                lenAlignment -= 1

    return perfMatch/lenAlignment * 100, missMatch/lenAlignment * 100, gapRef/lenAlignment * 100, gapOut/lenAlignment * 100


########################################################################################
########################################################################################
########################################################################################


orfology_tableHeader = ["OGID", "Ref_MSA_start", "Ref_MSA_end", "Ref_orf_start", "Ref_orf_end", "OutSpecies", "OutSpeciesORF", "OutS_MSA_start", "OutS_MSA_end", "OutS_orf_start", "OutS_orf_end", "NuclOverlap", "ClassOverlap"]

#ORF_table = pd.read_csv(orfology_table_path, names=orfology_tableHeader, sep="\t")
ORF_table = pd.read_csv(orfology_table_path, sep="\t")
ORF_table["refTranscript"] = ORF_table.RefORF.replace(":.*", "", regex = True)
ORF_table["refGene"] = ORF_table.refTranscript.replace("-T[1-9]", "", regex = True)
ORF_table["outTranscript"] = ORF_table.OutSpeciesORF.replace(":.*", "", regex = True)
ORF_table["outGene"] = ORF_table.outTranscript.replace("-T[1-9]", "", regex = True)
perIdentity = []
perRefGaps = []
perOutGaps = []

for index, row in ORF_table.iterrows():

    OGID = str(row["OGID"]).zfill(7) + ".msa.fa"
    MSA_path = MSA_dir_path + OGID
    sequences = SeqIO.index(MSA_path, "fasta")
    MSA_start = row["Ref_MSA_start"]
    MSA_end = row["Ref_MSA_end"]
    OutSp = row["OutSpecies"]
    sequenceRef = ""
    sequenceOutSp = ""
    for transcript in sequences:
        if sequenceRef == "":
            sequenceRef = str(sequences[transcript].seq).upper()
        elif transcript == row["outGene"]:
            sequenceOutSp = str(sequences[transcript].seq).upper()
        else:
            continue

    pedId, perMiss, perRefGap, perOutGap = GetPerIdPerGap(sequenceRef[MSA_start - 1:MSA_end - 1], sequenceOutSp[MSA_start - 1:MSA_end - 1])
    perIdentity.append(pedId)
    perRefGaps.append(perRefGap)
    perOutGaps.append(perOutGap)
    # Now check position by position if they are equal or not

ORF_table["perIdentity"] = perIdentity
ORF_table["perRefGap"] = perRefGaps
ORF_table["perOutGap"] = perOutGaps

ORF_table.to_csv(outPath, sep='\t', index=False)
