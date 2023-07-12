#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="This programs needs the config file to work with the desired species",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--ORFlist", default="", type=str, help="path to the file that has all the valid ORFs from the reference species, splitted by commas")
parser.add_argument("-OS", "--OutSpecies", default="", type=str, help="paths of the candidates files of the other species")
parser.add_argument("-orth", "--orthologs", default="", type=str, help="path to the tsv file that indicates the orthology of the genes")
parser.add_argument("-msa", "--MSAdir", default="", type=str, help="path to the directory where are all the regions aligned")
parser.add_argument("-o", "--outFile", default="", type=str, help="Name of the output file")

#parser.add_argument("-orf", "--ORF_type", default="", type=str, help="Provide information about the type of ORF that is going to be studied")
#parser.add_argument("-d", "--table_of_lengths", default="", type=str, help="Table that include the information about the length of each transcript and its parts")


#######################################################################################
############################### Inputs ################################################
#######################################################################################


args = parser.parse_args()

ScerORF_table_path = args.ORFlist
outOrfs_paths = args.OutSpecies.split(",")

Orthologues_table_path = args.orthologs
MSA_dir_path = args.MSAdir

if not MSA_dir_path.endswith("/"):
    MSA_dir_path = MSA_dir_path + "/"

ConservedORFs_path = args.outFile

# ScerORF_table_path = "/users/genomics/jmontanes/EvolutionaryNanopore/AdditionalSamplesRiboseq/riboNovel-Scer/Scer/Correct_format_files/candidateORF.genepred.fixed.noOv.txt"
# outOrfs_paths = "/users/genomics/jmontanes/EvolutionaryNanopore/AdditionalSamplesRiboseq/riboNovel-Spar/Spar/Correct_format_files/candidateORF.genepred.fixed.noOv.txt","/users/genomics/jmontanes/EvolutionaryNanopore/AdditionalSamplesRiboseq/riboNovel-Sbay/Sbay/Correct_format_files/candidateORF.genepred.fixed.noOv.txt"
#
# Orthologues_table_path = "/home/jmontanes/Documents/EvolutionNanopore/Outputs/OutputsR/One-to-one_orthologues.tsv"
# MSA_dir_path = "/home/jmontanes/Documents/EvolutionNanopore/Outputs/EvolutionNanopore/Outputs/Transcripts/transcriptAlignments/Fasta/"
# ConservedORFs_path = "Conserved_uORF.tsv"

genePredHeader = ["orfID", "chrom",	"strand", "start", "end", "ORFstart", "ORFend", "exonNumber", "StartExon", "EndExon"]
patternTranscriptInit = r':(\d+):'
patternTranscriptEnd = r'\|\d+:\d+:(\d+)\|'
patternTranscriptType = r'\:\d+\|(\w+)\|'
orfTypes = ["canonical","uORF", "ouORF", "dORF", "odORF"]

ORF_table = pd.read_csv(ScerORF_table_path, names=genePredHeader, sep="\t")
ORF_table["transcript"] = ORF_table.orfID.replace(":.*", "", regex = True)
ORF_table["gene"] = ORF_table.transcript.replace("-T[1-9]", "", regex = True)
ORF_table["transcriptInit"] = ORF_table.orfID.str.extract(patternTranscriptInit).astype(int)
ORF_table["transcriptEnd"] = ORF_table.orfID.str.extract(patternTranscriptEnd).astype(int) - 1
ORF_table["transcriptType"] = ORF_table.orfID.str.extract(patternTranscriptType)
ORF_table = ORF_table[ORF_table.transcriptType.isin(orfTypes)]

Orthologues_table = pd.read_csv(Orthologues_table_path, sep = "\t")
spFinder = {}
for sp in Orthologues_table.columns.to_list():
    spFinder[sp] = Orthologues_table[sp].to_list()

# get the list of pandas df

dfOut_list = list()
for file in outOrfs_paths:
    outORF_table = pd.read_csv(file, names=genePredHeader, sep="\t")
    outORF_table["transcript"] = outORF_table.orfID.replace(":.*", "", regex=True)
    outORF_table["gene"] = outORF_table.transcript.replace("-T[1-9]", "", regex=True)
    outORF_table["transcriptInit"] = outORF_table.orfID.str.extract(patternTranscriptInit).astype(int)
    outORF_table["transcriptEnd"] = outORF_table.orfID.str.extract(patternTranscriptEnd).astype(int) - 1
    outORF_table["transcriptType"] = outORF_table.orfID.str.extract(patternTranscriptType)
    outORF_table = outORF_table[outORF_table.transcriptType.isin(orfTypes)]

    spName = ""

    # Check which species is
    for transcript in outORF_table.transcript.replace("-T[1-9]", "", regex=True):
        for sp in spFinder.keys():
            if transcript in spFinder[sp]:
                spName = sp
                break
        if spName != "":
            break

    outORF_table["Species"] = spName
    dfOut_list.append(outORF_table)

########################################################################################
########################################################################################
########################################################################################

########################################################################################
################################ Functions #############################################
########################################################################################

### To extract positions in MSA

def GetMSAORFpos(nuclSeq,ORFStart,ORFEnd):
    MSApos = 0
    MSAstart = 0
    MSAend = 0
    ORFStart_local = ORFStart
    ORFEnd_local = ORFEnd

    for nucl in nuclSeq:
        if nucl == "-":
            MSApos += 1
        else:
            if ORFStart_local == 1:
                MSAstart = MSApos

            ORFStart_local = ORFStart_local - 1

            if ORFEnd_local == 1:
                MSAend = MSApos
                break

            ORFEnd_local = ORFEnd_local - 1
            MSApos += 1

    return MSAstart, MSAend

def GetORFpos(nuclSeq,MSAstart,MSAend):
    ORFpos = 0
    ORFStart_local = 0
    ORFEnd_local = 0

    for nucl in nuclSeq:
        if nucl == "-":

            if MSAstart == 1:
                ORFStart_local = ORFpos
            MSAstart = MSAstart - 1

            if MSAend == 1:
                ORFEnd_local = ORFpos
                break
            MSAend = MSAend - 1
        else:
            ORFpos += 1

            if MSAstart == 1:
                ORFStart_local = ORFpos
            MSAstart = MSAstart - 1
            if MSAend == 1:
                ORFEnd_local = ORFpos
                break
            MSAend = MSAend - 1


    return ORFStart_local, ORFEnd_local

def returnOverlap(RefMSAst, RefMSAend, SpMSAst, SpMSAend):
    if SpMSAst >= RefMSAst:
        if SpMSAst >= RefMSAst and SpMSAend < RefMSAend:
            # ORF inside the reference
            return ["Inner", SpMSAend - SpMSAst + 1]
        elif SpMSAend == RefMSAend and SpMSAst == RefMSAst:
            # ORF perfect match
            return ["Perfect match", RefMSAend - RefMSAst]
        else:
            # Right end larger
            return ["Right", RefMSAend - SpMSAst]
    else:
        if SpMSAend > RefMSAend:
            # Overlaps the ORF of the reference species
            return ["Outer", RefMSAend - RefMSAst]
        else:
            # Left and larger
            return ["Left", SpMSAend - RefMSAst]

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
    if lenAlignment > 0:
         return perfMatch/lenAlignment * 100, missMatch/lenAlignment * 100, gapRef/lenAlignment * 100, gapOut/lenAlignment * 100
    else:
        return 0, 100, 100, 100
########################################################################################
########################################################################################
########################################################################################

outOrthoORFs = 0
orthoORFs = 0

ConservedORFs = open(ConservedORFs_path, "w+")

# Lists

refORF = []
OGID_list = []
refMSA_startList = []
refMSA_endList = []
refTranscript_startList = []
refTranscript_endList = []

OutSpecies = []
OutSpeciesORF = []
OutSpeciesGene = []
OutSpeciesMSA_startList = []
OutSpeciesMSA_endList = []
OutSpeciesTranscript_startList = []
OutSpeciesTranscript_endList = []

ORFOverlap = []
ORFOverlapClass = []

perIdentity = []
perRefGaps = []
perOutGaps = []
for index, row in ORF_table.iterrows():

    orfID = row["orfID"]

    geneID, ribORFchr, strand_nORF_transcriptLength, refRibORFstart, ribORFend_ORFType_StartCodon = orfID.split(":")
    ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
    refRibORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")

    if ORFType not in ["canonical", "uORF", "ouORF", "dORF", "odORF"] :
        continue
    ####################
    # Manual fixing
    gene = geneID.split("-T")[0]
    OGID = Orthologues_table["OGID"][Orthologues_table["Scer"] == gene].to_string(index = False).zfill(7) + ".msa.fa"
    OGID_clean = OGID.split(".")[0]
    ####################
    ####################

    # Testing
    # if OGID_clean != "0004603":
    #    continue

    # The following means that it doesn't have orthologues
    if OGID.startswith("Series([]"):
        continue

    MSA_path = MSA_dir_path + OGID

    # There is a MSA alignment?
    try:
        sequences = SeqIO.index(MSA_path, "fasta")
    except:
        continue

    # The reference is always the first
    refMSAstart = 0
    refMSAend = 0

    sequenceRef = ""
    sequenceOutSp = ""

    for transcript in sequences:

        if refMSAstart == 0 and refMSAend == 0:
            refMSAstart, refMSAend = GetMSAORFpos(str(sequences[transcript].seq).upper(), int(refRibORFstart), int(refRibORFend))
            sequenceRef = str(sequences[transcript].seq).upper()
        else:
            for outDf in dfOut_list:
                outGenes = outDf.gene.unique().tolist()
                outTranscripts = outDf.transcript.unique().tolist()
                if transcript in outGenes:
                    ORFstart, ORFend = GetORFpos(str(sequences[transcript].seq).upper(), refMSAstart, refMSAend)

                    if len(list(set(str(sequences[transcript].seq).upper()[refMSAstart:refMSAend]))) == 1:
                        continue

                    # Filter main df

                    mask = ((outDf.transcriptEnd > ORFstart) & (outDf.transcriptInit < ORFend) & (outDf.gene == transcript))
                    subDf = outDf[mask].copy()

                    if len(subDf) == 0:
                        continue


                    for i in subDf.index.to_list():
                        # First add the information of the reference

                        outSpGeneID, outSpRibORFchr, outSpStrand_nORF_transcriptLength, outSpRibORFstart, outSpRibORFend_ORFType_StartCodon = subDf["orfID"][i].split(":")
                        outSpRibORFstrand, outSpnORF, outSpTranscriptLength = outSpStrand_nORF_transcriptLength.split("|")
                        outSpRibORFend, outSpORFType, outSpStartCodon = outSpRibORFend_ORFType_StartCodon.split("|")

                        outSpMSAstart, outSpMSAend = GetMSAORFpos(str(sequences[transcript].seq).upper(), int(outSpRibORFstart) , int(outSpRibORFend) - 1)

                        ClassOv, NuclOv = returnOverlap(refMSAstart, refMSAend, outSpMSAstart, outSpMSAend + 1)

                        refORF.append(orfID)
                        OGID_list.append(OGID_clean)
                        refMSA_startList.append(refMSAstart + 1)
                        refMSA_endList.append(refMSAend + 1)
                        refTranscript_startList.append(int(refRibORFstart))
                        refTranscript_endList.append(int(refRibORFend) - 1)

                        OutSpeciesMSA_startList.append(outSpMSAstart + 1)
                        OutSpeciesMSA_endList.append(outSpMSAend + 2)
                        OutSpeciesTranscript_startList.append(outSpRibORFstart)
                        OutSpeciesTranscript_endList.append(int(outSpRibORFend) - 1)

                        # Information about the out Species overlap

                        OutSpecies.append(subDf["Species"][i])
                        OutSpeciesORF.append(subDf["orfID"][i])

                        # if subDf["Species"][i] == "Sbay":
                        #     sbaySeq = str(sequences[transcript].seq).upper()

                        ORFOverlap.append(NuclOv)
                        ORFOverlapClass.append(ClassOv)

                        OutSpeciesGene.append(transcript)
                        sequenceOutSp = str(sequences[transcript].seq).upper()
                        pedId, perMiss, perRefGap, perOutGap = GetPerIdPerGap(sequenceRef[refMSAstart:refMSAend], sequenceOutSp[refMSAstart:refMSAend])
                        perIdentity.append(pedId)
                        perRefGaps.append(perRefGap)
                        perOutGaps.append(perOutGap)

                        # print(subDf["orfID"][i])
                        # print(ORFstart)
                        # print(ORFend)
                        # print(refMSAstart)
                        # print(refMSAend)




# refSeq = sequenceRef

# Create a set of all unique strings from both lists

df = pd.DataFrame(list(zip(refORF, OGID_list, refMSA_startList, refMSA_endList, refTranscript_startList, refTranscript_endList,
                           OutSpecies, OutSpeciesORF, OutSpeciesMSA_startList, OutSpeciesMSA_endList, OutSpeciesTranscript_startList, OutSpeciesTranscript_endList,
                           ORFOverlap, ORFOverlapClass,
                           perIdentity, perRefGaps, perOutGaps
                           )), columns=['RefORF', 'OGID', 'Ref_MSA_start', 'Ref_MSA_end', 'Ref_orf_start', 'Ref_orf_end',
                                        'OutSpecies', 'OutSpeciesORF', 'OutS_MSA_start', 'OutS_MSA_end', 'OutS_orf_start', 'OutS_orf_end',
                                        "NuclOverlap","ClassOverlap",
                                        "perIdentity", "perRefGap", "perOutGap"])
# Output the DataFrame to TSV format
df.to_csv(ConservedORFs_path, sep='\t', index=False)
