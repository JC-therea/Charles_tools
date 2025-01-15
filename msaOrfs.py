#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import argparse

def readOrthologyTable(orthologyTablePath):
    Orthologues_table = pd.read_csv(Orthologues_table_path, sep = "\t")
    spFinder = {}
    for sp in Orthologues_table.columns.to_list():
        spFinder[sp] = Orthologues_table[sp].to_list()
    return Orthologues_table, spFinder

def readGenePred(gPredPath, orthologyDictObject):

    genePredHeader = ["orfID", "chrom",	"strand", "start", "end", "ORFstart", "ORFend", "exonNumber", "StartExon", "EndExon"]
    patternTranscriptInit = r':(\d+):'
    patternTranscriptEnd = r'\|\d+:\d+:(\d+)\|'
    patternTranscriptType = r'\:\d+\|(\w+)\|'
    orfTypes = ["canonical","uORF", "ouORF", "dORF", "odORF"]

    ORF_table = pd.read_csv(gPredPath, names=genePredHeader, sep="\t")
    ORF_table["transcript"] = ORF_table.orfID.replace(":.*", "", regex = True)
    ORF_table["transcriptInit"] = ORF_table.orfID.str.extract(patternTranscriptInit).astype(int)
    ORF_table["transcriptEnd"] = ORF_table.orfID.str.extract(patternTranscriptEnd).astype(int)
    ORF_table["transcriptType"] = ORF_table.orfID.str.extract(patternTranscriptType)
    ORF_table = ORF_table[ORF_table.transcriptType.isin(orfTypes)]

    spName = ""

    # Check which species is
    for transcript in ORF_table["transcript"]:
        if spName != "":
            break
        for sp in orthologyDictObject.keys():
            if transcript in orthologyDictObject[sp]:
                spName = sp
                break

    ORF_table["Species"] = spName

    return ORF_table

def coordTranscript2Msa(nuclSeq,ORFStart,ORFEnd):
    MSApos = 0
    MSAstart = 0
    MSAend = 0
    ORFStart_local = int(ORFStart)
    ORFEnd_local = int(ORFEnd)

    for nucl in nuclSeq:
        if nucl == "-":
            MSApos += 1
        else:
            if ORFStart_local == 1:
                MSAstart = MSApos + 1

            ORFStart_local = ORFStart_local - 1

            if ORFEnd_local == 1:
                MSAend = MSApos
                break

            ORFEnd_local = ORFEnd_local - 1
            MSApos += 1
    if MSAend == 0:
        MSAend = MSApos
    return MSAstart, MSAend

def coordMsa2Transcript(nuclSeq,MSAstart,MSAend):
    ORFpos = 0
    ORFStart_local = 0
    ORFEnd_local = 0
    MSAstart = int(MSAstart)
    MSAend = int(MSAend)

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
        if SpMSAend == RefMSAend and SpMSAst == RefMSAst:
            return ["Perfect match", RefMSAend - RefMSAst +1 ]
        elif SpMSAend <= RefMSAend:
            # ORF inside the reference
            return ["Inner", SpMSAend - SpMSAst +1]
        elif SpMSAend == RefMSAend and SpMSAst == RefMSAst:
            # ORF perfect match
            return ["Perfect match", RefMSAend - RefMSAst +1]
        else:
            # Right end larger
            return ["Right", RefMSAend - SpMSAst +1]
    else:
        if SpMSAend > RefMSAend:
            # Overlaps the ORF of the reference species
            return ["Outer", RefMSAend - RefMSAst +1]
        else:
            # Left and larger
            return ["Left", SpMSAend - RefMSAst +1 ]

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

def msaOrfs(genePredsPaths,orthologTablePath,MsaDirPath,outFilePath):

    with open(outFilePath, "w+") as outfile:

        orthologyTable, speciesDict = readOrthologyTable(orthologTablePath)
        refGenePredPath = genePredsPaths.split(",")[0]
        refGenePred = readGenePred(refGenePredPath, speciesDict)

        outGenePredPaths = genePredsPaths.split(",")[1:]
        outGenePredList = []

        for outGenePredPath in outGenePredPaths:
            gPred = readGenePred(outGenePredPath, speciesDict)
            outGenePredList.append(gPred)

        # The order is supposed to be the same that the one in the MSA
        # Print headers
        header = ['RefORF', 'OGID', 'Ref_MSA_start', 'Ref_MSA_end', 'Ref_orf_start', 'Ref_orf_end',
                  'OutSpecies', 'OutSpeciesORF', 'OutS_MSA_start', 'OutS_MSA_end', 'OutS_orf_start', 'OutS_orf_end',
                  "NuclOverlap","ClassOverlap", "perIdentity", "perRefGap", "perOutGap"]

        outfile.write("\t".join(map(str, header)) + "\n")
        # Now iterate over each orthologue in orthologyTable

        nOrthologs = len(orthologyTable)
        perToShow = 5
        rowCounter = 0.0
        for index, row in orthologyTable.iterrows():
            rowCounter += 1
            if rowCounter/nOrthologs * 100 > perToShow:
                print(f'More than {perToShow}% of sequences analyzed')
                perToShow += 5
            OGID = str(row['OGID']).zfill(7)
            msaPath = MsaDirPath + OGID + ".msa.fa"
            sequences = SeqIO.index(msaPath, "fasta")

            # Get all the orfs in the reference sequence

            refGenePredSubset = refGenePred.loc[(refGenePred['transcript']) == row.iloc[0]]
            outGenePredListSubset = []

            # Get all the orfs in the outgroupspecies

            for i in range(1, len(row) - 1):

                outPos = i - 1
                outGPred = outGenePredList[outPos]
                outTranscript = row.iloc[i]

                # Exception for some coordTranscripts
                if outTranscript.startswith("rna-XM"):
                    outTranscript = outTranscript.split("rna-")[1]
                outGPredSubset = outGPred[ outGPred['transcript'] == outTranscript]
                outGenePredListSubset.append(outGPredSubset)

            # Now check orf by orf of the reference if there is some overlap
            for indexOrf, rowOrf in refGenePredSubset.iterrows():
                    orfID = rowOrf["orfID"]

                    transcriptId, ribORFchr, strand_nORF_transcriptLength, refRibORFstart, ribORFend_ORFType_StartCodon = orfID.split(":")
                    ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
                    refRibORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")

                    refRibORFend = int(refRibORFend)

                    sequenceIds = list(sequences.keys())
                    firstSequenceId = sequenceIds[0]
                    refSeq = sequences[firstSequenceId].seq.upper()

                    refStartMsa, refEndMsa = coordTranscript2Msa(refSeq, refRibORFstart, refRibORFend)

                    for outI in range(1, len(sequenceIds) ):

                        outPos = outI - 1
                        outSequenceId = sequenceIds[outI]

                        if outSequenceId == "*":
                            continue

                        outSeq = sequences[outSequenceId].seq.upper()
                        outGPred = outGenePredListSubset[outPos]
                        outSpecies = list(speciesDict.keys())[outI]

                        coordsOutStartOrf, coordsOutEndOrf = coordMsa2Transcript(outSeq, refStartMsa, refEndMsa)

                        innerOutGPred = outGPred[((outGPred["transcriptInit"] >= coordsOutStartOrf) & (outGPred["transcriptInit"] <= coordsOutEndOrf) ) |
                                                 ((outGPred["transcriptEnd"] > coordsOutStartOrf) & (outGPred["transcriptEnd"] < coordsOutEndOrf))
                                            ]
                        if innerOutGPred.empty:
                            continue
                        for indexOutOrf, rowOutOrf in innerOutGPred.iterrows():
                            # First add the information of the reference
                            outOrf = rowOutOrf["orfID"]
                            outStartOrf = int(rowOutOrf["transcriptInit"])
                            outEndOrf = int(rowOutOrf["transcriptEnd"])
                            outStartMsa, outEndMsa = coordTranscript2Msa(outSeq, outStartOrf, outEndOrf)
                            overlapType, nuclOverlapping = returnOverlap(refStartMsa, refEndMsa, outStartMsa, outEndMsa)
                            pedId, perMiss, perRefGap, perOutGap = GetPerIdPerGap(refSeq[refStartMsa:refEndMsa], outSeq[refStartMsa:refEndMsa])


                            outfile.write(f"{orfID}\t{OGID}\t{refStartMsa}\t{refEndMsa}\t{refRibORFstart}\t{refRibORFend}\t{outSpecies}\t{outOrf}\t{outStartMsa}\t{outEndMsa}\t{outStartOrf}\t{outEndOrf}\t{nuclOverlapping}\t{overlapType}\t{pedId}\t{perRefGap}\t{perOutGap}\n")


# Program starts here


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This program is used to extract orfs aligned",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--ORFlists", default="", type=str, help="path of all the genePredFiles, first will be used as reference")
    parser.add_argument("-orth", "--orthologs", default="", type=str, help="path to the tsv file that indicates the orthology of the genes")
    parser.add_argument("-msa", "--MSAdir", default="", type=str, help="path to the directory where are all the regions aligned")
    parser.add_argument("-o", "--outFile", default="", type=str, help="Name of the output file")

    #######################################################################################
    ############################### Inputs ################################################
    #######################################################################################


    args = parser.parse_args()
    genePredPaths = args.ORFlists
    Orthologues_table_path = args.orthologs
    MSA_dir_path = args.MSAdir

    if not MSA_dir_path.endswith("/"):
        MSA_dir_path = MSA_dir_path + "/"

    ConservedORFs_path = args.outFile

    msaOrfs(genePredPaths,Orthologues_table_path,MSA_dir_path,ConservedORFs_path)
