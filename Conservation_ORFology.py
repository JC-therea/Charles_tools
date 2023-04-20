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

for index, row in ORF_table.iterrows():

    orfID = row["orfID"]
    gene = row["gene"]

    #spsAlign = Orthologues_table["speciesAligment"][Orthologues_table["Scer"] == gene].to_string(index = False)
    OGID = Orthologues_table["OGID"][Orthologues_table["Scer"] == gene].to_string(index = False).zfill(7) + ".mfa"
    if OGID.startswith("Series([]"):
        continue
    MSA_path = MSA_dir_path + OGID

    geneID, ribORFchr, strand_nORF_transcriptLength, ribORFstart, ribORFend_ORFType_StartCodon = orfID.split(":")
    ribORFstrand, nORF, transcriptLength = strand_nORF_transcriptLength.split("|")
    ribORFend, ORFType, StartCodon = ribORFend_ORFType_StartCodon.split("|")

    transcriptID_list = []
    try:
        sequences = SeqIO.index(MSA_path, "fasta")
    except:
        continue

    # S.cerevisiae is always the first
    MSAstart = 0
    MSAend = 0
    Spar = []
    Spar_gene = ""
    Sbay = []
    Sbay_gene = ""
    for transcript in sequences:
        #print(transcript)
        if MSAstart == 0 and MSAend == 0:
            MSAstart, MSAend = GetMSAORFpos(str(sequences[transcript].seq).upper(), int(ribORFstart), int(ribORFend))
        elif transcript in Spar_genes:
            ORFstart, ORFend = GetORFpos(str(sequences[transcript].seq).upper(), MSAstart, MSAend)
            Spar.append(ORFstart)
            Spar.append(ORFend)
            Spar_gene = transcript
        elif transcript in Sbay_genes:
            ORFstart, ORFend = GetORFpos(str(sequences[transcript].seq).upper(), MSAstart, MSAend)
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

        if len(SparORF_filtered_df) >= 1:
            Spar_orthoORFs += 1
            ORF_in_Spar.append(orfID)
            Spar_ORFs[orfID] = "{-}".join(SparORF_filtered_df["orfID"])

    if len(Sbay_gene) > 2 and sum(SbayORF_table.gene == Sbay_gene) > 0:
        SbayORF_table_subset = SbayORF_table[SbayORF_table.gene == Sbay_gene].copy()
        SbayORF_table_subset['within_window'] = (SbayORF_table_subset['transcriptInit'] >= Sbay[0]) & (
                SparORF_table_subset['transcriptEnd'] <= Sbay[1])

        # Create a new column to store whether each transcript is partially within the window
        SbayORF_table_subset['partially_within_window'] = ((SbayORF_table_subset['transcriptInit'] < Sbay[0]) & (
                SbayORF_table_subset['transcriptEnd'] > Sbay[0])) | (
                                                                  (SbayORF_table_subset['transcriptInit'] < Sbay[
                                                                      1]) & (SbayORF_table_subset['transcriptEnd'] >
                                                                             Sbay[1]))

        # Filter the transcripts that are either completely within the window or partially within the window with at least 50% of their length inside the window
        SbayORF_filtered_df = SbayORF_table_subset[
            ((SbayORF_table_subset['within_window']) | (SbayORF_table_subset['partially_within_window']))].copy()

        if len(SbayORF_filtered_df) >= 1:
            Sbay_orthoORFs += 1
            ORF_in_Sbay.append(orfID)
            Sbay_ORFs[orfID] = "{-}".join(SbayORF_filtered_df["orfID"])

# Create a set of all unique strings from both lists
unique_transcripts = set(ORF_in_Spar + ORF_in_Sbay)
Spar_listTranscripts = []

Sbay_listTranscripts = []

for transcript in unique_transcripts:

    if transcript in Spar_ORFs.keys():
        Spar_listTranscripts.append(Spar_ORFs[transcript])
    else:
        Spar_listTranscripts.append("-")


    if transcript in Sbay_ORFs.keys():
        Sbay_listTranscripts.append(Sbay_ORFs[transcript])
    else:
        Sbay_listTranscripts.append("-")

df = pd.DataFrame(list(zip(unique_transcripts, Spar_listTranscripts, Sbay_listTranscripts)), columns=['ORF', 'Spar_ORFs', 'Sbay_ORFs'])

# Output the DataFrame to TSV format
df.to_csv(ConserveduORFs_path, sep='\t', index=False)
