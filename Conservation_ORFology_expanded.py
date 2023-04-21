import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="This programs needs the config file to work with the desired species",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--ORFlist", default="", type=str, help="path to the file that has all the valid ORFs from the reference species")
parser.add_argument("-spar", "--ORFSpar", default="", type=str, help="path to the file that has all the ORFs from the alt species")
parser.add_argument("-sbay", "--ORFSbay", default="", type=str, help="path to the file that has all the ORFs from the alt species")
parser.add_argument("-orth", "--orthologs", default="", type=str, help="path to the tsv file that indicates the orthology of the genes")
parser.add_argument("-msa", "--MSAdir", default="", type=str, help="path to the directory where are all the regions aligned")
parser.add_argument("-o", "--outFile", default="", type=str, help="Name of the output file")

#parser.add_argument("-orf", "--ORF_type", default="", type=str, help="Provide information about the type of ORF that is going to be studied")
#parser.add_argument("-d", "--table_of_lengths", default="", type=str, help="Table that include the information about the length of each transcript and its parts")

args = parser.parse_args()

ScerORF_table_path = args.ORFlist
SparORF_table_path = args.ORFSpar
SbayORF_table_path = args.ORFSbay

Orthologues_table_path = args.orthologs
MSA_dir_path = args.MSAdir

if not MSA_dir_path.endswith("/"):
    MSA_dir_path = MSA_dir_path + "/"

ConserveduORFs_path = args.outFile

# ScerORF_table_path = "/home/jmontanes/Documents/EvolutionNanopore/Outputs/OutputsR/RibosomeProfilingStats/Scer_validuORFs.tsv"
# SparORF_table_path = "/home/jmontanes/Documents/EvolutionNanopore/Outputs/OutputsR/RibosomeProfilingStats/Spar_ORFs.tsv"
# SbayORF_table_path = "/home/jmontanes/Documents/EvolutionNanopore/Outputs/OutputsR/RibosomeProfilingStats/Sbay_ORFs.tsv"

# Orthologues_table_path = "/home/jmontanes/Documents/EvolutionNanopore/Outputs/OutputsR/Evolution/ORFevolution/Orthology_evolution_5utr.tsv"
# MSA_dir_path = "/home/jmontanes/Documents/EvolutionNanopore/Outputs/EvolutionNanopore/Outputs/Transcripts/UTR5/UTR5_alignments/"
# ConserveduORFs_path = "Conserved_uORF.txt"

ORF_table = pd.read_csv(ScerORF_table_path, sep = "\t")
SparORF_table = pd.read_csv(SparORF_table_path, sep = "\t")
SbayORF_table = pd.read_csv(SbayORF_table_path, sep = "\t")

Orthologues_table = pd.read_csv(Orthologues_table_path, sep = "\t")
Spar_genes = Orthologues_table["Spar"].to_list()
Sbay_genes = Orthologues_table["Sbay"].to_list()

# Test
### Fun to extract positions in MSA

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
            MSApos += 1
            ORFStart_local = ORFStart_local - 1
            ORFEnd_local = ORFEnd_local - 1

            if ORFStart_local == 1:
                MSAstart = MSApos
            if ORFEnd_local == 1:
                MSAend = MSApos
                break

    return MSAstart, MSAend

def GetORFpos(nuclSeq,MSAstart,MSAend):
    ORFpos = 0
    ORFStart_local = 0
    ORFEnd_local = 0

    for nucl in nuclSeq:
        if nucl == "-":
            MSAstart = MSAstart - 1
            MSAend = MSAend - 1

            if MSAstart == 1:
                ORFStart_local = ORFpos
            if MSAend == 1:
                ORFEnd_local = ORFpos
                break
        else:
            ORFpos += 1
            MSAstart = MSAstart - 1
            MSAend = MSAend - 1
            if MSAstart == 1:
                ORFStart_local = ORFpos
            if MSAend == 1:
                ORFEnd_local = ORFpos
                break

    return ORFStart_local, ORFEnd_local

Spar_orthoORFs = 0
Sbay_orthoORFs = 0
orthoORFs = 0

ConserveduORFs = open(ConserveduORFs_path, "w+")
ORF_in_Spar = []
Spar_ORFs = {}
ORF_in_Sbay = []
Sbay_ORFs = {}

# Lists
ScerOrfIDList = []
OGID_list = []
ScerMSA_start_list = []
ScerMSA_end_list = []
ScerStartORFinT_list = []
ScerEndORFinT_list = []
OutSpecies = []
OutSpeciesORFid = []
OutSpeciesMSAStart_list = []
OutSpeciesMSAEnd_list = []
OutSpeciesORFstart_list = []
OutSpeciesORFEnd_list = []

for index, row in ORF_table.iterrows():

    orfID = row["orfID"]
    gene = row["gene"]
    #spsAlign = Orthologues_table["speciesAligment"][Orthologues_table["Scer"] == gene].to_string(index = False)
    OGID = Orthologues_table["OGID"][Orthologues_table["Scer"] == gene].to_string(index = False).zfill(7) + ".mfa"
    OGID_clean = OGID.split(".")[0]
    if OGID.startswith("Series([]"):
        continue
    MSA_path = MSA_dir_path + OGID

    geneID, ribORFchr, strand_nORF_transcriptLength, ScerribORFstart, ribORFend_ORFType_StartCodon = orfID.split(":")
    ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
    ScerribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")

    transcriptID_list = []
    try:
        sequences = SeqIO.index(MSA_path, "fasta")
    except:
        continue

    # S.cerevisiae is always the first
    ScerMSAstart = 0
    ScerMSAend = 0
    Spar = []
    Spar_gene = ""
    Sbay = []
    Sbay_gene = ""
    for transcript in sequences:
        #print(transcript)
        if ScerMSAstart == 0 and ScerMSAend == 0:
            ScerMSAstart, ScerMSAend = GetMSAORFpos(str(sequences[transcript].seq).upper(), int(ScerribORFstart), int(ScerribORFend))
        elif transcript in Spar_genes:
            ORFstart, ORFend = GetORFpos(str(sequences[transcript].seq).upper(), ScerMSAstart, ScerMSAend)
            Spar.append(ORFstart)
            Spar.append(ORFend)
            Spar_gene = transcript
        elif transcript in Sbay_genes:
            ORFstart, ORFend = GetORFpos(str(sequences[transcript].seq).upper(), ScerMSAstart, ScerMSAend)
            Sbay.append(ORFstart)
            Sbay.append(ORFend)
            Sbay_gene = transcript

    # Now check if there is ORFs from both species
    # S.par
    if len(Spar_gene) > 2 and sum(SparORF_table.gene == Spar_gene) > 0:
        SparORF_table_subset = SparORF_table[SparORF_table.gene == Spar_gene].copy()
        SparORF_table_subset['within_window'] = (SparORF_table_subset['transcriptInit'] >= Spar[0]) & (SparORF_table_subset['transcriptEnd'] <= Spar[1])

        # Create a new column to store whether each transcript is partially within the window
        SparORF_table_subset['partially_within_window'] = ((SparORF_table_subset['transcriptInit'] < Spar[0]) & (SparORF_table_subset['transcriptEnd'] > Spar[0])) | (
                    (SparORF_table_subset['transcriptInit'] < Spar[1]) & (SparORF_table_subset['transcriptEnd'] > Spar[1]))

        # Filter the transcripts that are either completely within the window or partially within the window with at least 50% of their length inside the window
        SparORF_filtered_df = SparORF_table_subset[((SparORF_table_subset['within_window']) | (SparORF_table_subset['partially_within_window']))].copy()
        SparORF_filtered_df = SparORF_filtered_df.reset_index()

        for i in range(0,len(SparORF_filtered_df)):
            ScerOrfIDList.append(orfID)
            OGID_list.append(OGID_clean)
            ScerMSA_start_list.append(ScerMSAstart + 1)
            ScerMSA_end_list.append(ScerMSAend + 1)
            ScerStartORFinT_list.append(ScerribORFstart)
            ScerEndORFinT_list.append(ScerribORFend)
            OutSpecies.append("S.paradoxus")

        if len(SparORF_filtered_df) >= 1:
            Spar_orthoORFs += 1
            ORF_in_Spar.append(orfID)
            Spar_ORFs[orfID] = "{-}".join(SparORF_filtered_df["orfID"])
            for i in range(0, len(SparORF_filtered_df)):
                OutSpeciesORFid.append(SparORF_filtered_df["orfID"][i])
                geneID, ribORFchr, strand_nORF_transcriptLength, ribORFstart, ribORFend_ORFType_StartCodon = SparORF_filtered_df["orfID"][i].split(
                    ":")
                ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
                ribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")
                # From here extract the MSA position
                for transcript in sequences:
                    if transcript in Spar_genes:
                        MSAstart, MSAend = GetMSAORFpos(str(sequences[transcript].seq).upper(), int(ribORFstart), int(ribORFend))
                OutSpeciesMSAStart_list.append(MSAstart + 1)
                OutSpeciesMSAEnd_list.append(MSAend + 1)
                OutSpeciesORFstart_list.append(ribORFstart)
                OutSpeciesORFEnd_list.append(ribORFend)
    # Sbay
    if len(Sbay_gene) > 2 and sum(SbayORF_table.gene == Sbay_gene) > 0:
        SbayORF_table_subset = SbayORF_table[SbayORF_table.gene == Sbay_gene].copy()
        SbayORF_table_subset['within_window'] = (SbayORF_table_subset['transcriptInit'] >= Sbay[0]) & (
                SbayORF_table_subset['transcriptEnd'] <= Sbay[1])

        # Create a new column to store whether each transcript is partially within the window
        SbayORF_table_subset['partially_within_window'] = ((SbayORF_table_subset['transcriptInit'] < Sbay[0]) & (
                SbayORF_table_subset['transcriptEnd'] > Sbay[0])) | (
                                                                  (SbayORF_table_subset['transcriptInit'] < Sbay[
                                                                      1]) & (SbayORF_table_subset['transcriptEnd'] >
                                                                             Sbay[1]))

        # Filter the transcripts that are either completely within the window or partially within the window with at least 50% of their length inside the window
        SbayORF_filtered_df = SbayORF_table_subset[
            ((SbayORF_table_subset['within_window']) | (SbayORF_table_subset['partially_within_window']))].copy()
        SbayORF_filtered_df = SbayORF_filtered_df.reset_index()
        for i in range(0,len(SbayORF_filtered_df)):
            ScerOrfIDList.append(orfID)
            OGID_list.append(OGID_clean)
            ScerMSA_start_list.append(ScerMSAstart + 1)
            ScerMSA_end_list.append(ScerMSAend + 1)
            ScerStartORFinT_list.append(ScerribORFstart)
            ScerEndORFinT_list.append(ScerribORFend)
            OutSpecies.append("S.bayanus")

        if len(SbayORF_filtered_df) >= 1:

            Sbay_orthoORFs += 1
            ORF_in_Sbay.append(orfID)
            Sbay_ORFs[orfID] = "{-}".join(SbayORF_filtered_df["orfID"])
            for i in range(0, len(SbayORF_filtered_df)):
                OutSpeciesORFid.append(SbayORF_filtered_df["orfID"][i])

                geneID, ribORFchr, strand_nORF_transcriptLength, ribORFstart, ribORFend_ORFType_StartCodon = SbayORF_filtered_df["orfID"][i].split(
                    ":")
                ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
                ribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")
                # From here extract the MSA position
                for transcript in sequences:
                    if transcript in Sbay_genes:
                        MSAstart, MSAend = GetMSAORFpos(str(sequences[transcript].seq).upper(), int(ribORFstart), int(ribORFend))
                OutSpeciesMSAStart_list.append(MSAstart + 1)
                OutSpeciesMSAEnd_list.append(MSAend + 1)
                OutSpeciesORFstart_list.append(ribORFstart)
                OutSpeciesORFEnd_list.append(ribORFend)

# Create a set of all unique strings from both lists
df = pd.DataFrame(list(zip(ScerOrfIDList, OGID_list, ScerMSA_start_list, ScerMSA_end_list, ScerStartORFinT_list, ScerEndORFinT_list,
                           OutSpecies,OutSpeciesORFid, OutSpeciesMSAStart_list, OutSpeciesMSAEnd_list, OutSpeciesORFstart_list, OutSpeciesORFEnd_list
                           )), columns=['RefORF', 'OGID', 'Ref_MSA_start', 'Ref_MSA_end', 'Ref_orf_start', 'Ref_orf_end',
                                        'OutSpecies', 'OutSpeciesORF', 'OutS_MSA_start', 'OutS_MSA_end', 'OutS_orf_start', 'OutS_orf_end'])

# Output the DataFrame to TSV format
df.to_csv(ConserveduORFs_path, sep='\t', index=False)
