
# FULL NEW SCRIPT TO EXTRACT ISOFORMS
import pandas as pd
import csv
import re
import argparse

# Standard inputs

parser = argparse.ArgumentParser(description="This programs extract the isoforms from the output of stringtie",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-st", "--StringTieAnnot", default="0", type=str, help="path to the annotation file created by Stringtie")
parser.add_argument("-d", "--defaultAnnotation", default="0", type=str, help="path to the annotation file created that was used with stringtie in gff")
parser.add_argument("-o", "--outDir", default="0", type=str, help="Output path file created by Stringtie")

args = parser.parse_args()

STGTFpath = args.StringTieAnnot
defaultAnnotPath = args.defaultAnnotation
OutDir = args.outDir

if not OutDir.endswith("/"):
    OutDir = OutDir + "/"

AllIsoPath = OutDir + "Any_isoforms.gtf"
outCDSpath = OutDir + "Any_isoforms_CDS.gtf"
IsoTranscriptsPath = OutDir + "New_isoforms_transcripts.gtf"
IsoCDSsPath = OutDir + "New_isoforms_CDS.gtf"

####### Custom inputs

# STGTFpath = "/home/jmontanes/Documents/Nanopore/Outputs/Hybrid_reads/GTF/Fmlrc_TranscriptClean_MAPQ_correction/Spom_stringtie_correction_ONT_with_annotation.gtf"
# defaultAnnotPath = "/home/jmontanes/Documents/Nanopore/Inputs/PomBase/Schizosaccharomyces_pombe_all_chromosomes.gff3"
# AllIsoPath = "/home/jmontanes/Documents/Software/Charles_tools/Charles_tools/Test/Any_isoforms.gtf"
# outCDSpath = "/home/jmontanes/Documents/Software/Charles_tools/Charles_tools/Test/Any_isoforms_CDS.gtf"
#
# IsoTranscriptsPath = "/home/jmontanes/Documents/Software/Charles_tools/Charles_tools/Test/New_isoforms_transcripts.gtf"
# IsoCDSsPath = "/home/jmontanes/Documents/Software/Charles_tools/Charles_tools/Test/New_isoforms_CDS.gtf"

#### First we keep those genes with at least one isoform.

STGTF = pd.read_table(STGTFpath, sep='\t', comment="#", dtype='str', header=None)
STGTF.columns = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
Gene_name = STGTF["attributes"].str.replace('gene_id "', '', regex=True)
Gene_name = Gene_name.str.replace('"; transcript_id.*', '', regex=True)
STGTF["Gene"] = Gene_name

STGTF_just_transcripts = STGTF[:][STGTF.type.isin(["transcript"])]
Genes = STGTF_just_transcripts["Gene"].tolist()
Transcripts_per_gene = {i: Genes.count(i) for i in Genes}
Genes_with_isoforms = []
for key, value in Transcripts_per_gene.items():
    if value > 1:
        Genes_with_isoforms.append(key)

STGTF_onlyIso = STGTF[:][STGTF.Gene.isin(Genes_with_isoforms)].copy()
STGTF_onlyIsoOut = STGTF_onlyIso.copy()
STGTF_onlyIsoOut = STGTF_onlyIsoOut.drop(columns = "Gene")
pd.DataFrame(STGTF_onlyIsoOut).to_csv(AllIsoPath, sep = "\t", index = False, header = False, quoting= csv.QUOTE_NONE)

#### Now we filter those isoforms that only modify the UTRs but not the CDS

Transcript_name = STGTF_onlyIso["attributes"].str.replace('.*transcript_id "', '', regex=True)
Transcript_name = Transcript_name.str.replace('"; .*', '', regex=True)

STGTF_onlyIso["Transcript"] = Transcript_name
STGTF_onlyIso["Ref_gene"] = "-"

Gene_name_u = STGTF_onlyIso["Gene"].unique().tolist()

for index, row in STGTF_onlyIso.iterrows():
    ref_ID = len(re.findall("ref_gene_id",row["attributes"]))
    if ref_ID > 0:
        ref_name = re.sub('.*ref_gene_id "', '', row["attributes"])
        ref_name = re.sub('"; .*', '', ref_name)
        row["Ref_gene"] = ref_name

# Now load the default annotation and filter by CDS

################################################################################################
# If you want to introduce a gtf file as reference please modify the following

defaultAnnot = pd.read_table(defaultAnnotPath, sep='\t', comment="#", dtype='str', header=None)
defaultAnnot.columns = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
defaultAnnotCDS = defaultAnnot[defaultAnnot.type == "CDS"].copy()
Gene_name = defaultAnnotCDS["attributes"].str.replace('ID=', '', regex=True)
Gene_name = Gene_name.str.replace('\..:.*', '', regex=True)
defaultAnnotCDS["Gene"] = Gene_name

# Record ncRNA
ncRNA_Genes = defaultAnnot["attributes"][defaultAnnot.type.isin(["ncRNA"])].str.replace('.*Parent=', '', regex=True).str.replace('\..:.*', '', regex=True)

################################################################################################
CDS_isoforms_list_of_lists = []

for gene in Gene_name_u:
    STGTF_onlyIso_subset = STGTF_onlyIso[:][STGTF_onlyIso.Gene.isin([gene])]
    Ref_gene_subset = STGTF_onlyIso_subset["Ref_gene"].unique().tolist()
    if "-" in Ref_gene_subset:
        Ref_gene_subset.remove("-")
    if len(Ref_gene_subset) == 0 or len(Ref_gene_subset) > 1:
        continue
    Ref_gene_subset_str = Ref_gene_subset[0]

    defaultAnnotCDS_subset = defaultAnnotCDS[:][defaultAnnotCDS.Gene.isin([Ref_gene_subset_str])]
    CDS_coordinates = defaultAnnotCDS_subset["start"].unique().tolist() + defaultAnnotCDS_subset["end"].unique().tolist()
    CDS_coordinates = [int(item) for item in CDS_coordinates]
    CDS_coordinates.sort()
    start_CDS = CDS_coordinates[0]
    end_CDS = CDS_coordinates[-1]

    Transcripts = STGTF_onlyIso_subset["Transcript"].unique().tolist()

    for transcript in Transcripts:
        STGTF_onlyIso_sub_subset = STGTF_onlyIso_subset[:][STGTF_onlyIso_subset.Transcript.isin([transcript])]
        for index, row in STGTF_onlyIso_sub_subset.iterrows():
            if row[2] != "exon":
                continue
            if int(row["end"]) < start_CDS:
                continue
            if int(row["start"]) > end_CDS:
                continue
            if int(row["start"]) < start_CDS:
                row["start"] = start_CDS
            if int(row["end"]) > end_CDS:
                row["end"] = end_CDS
            row["type"] = "CDS"
            CDS_isoforms_list_of_lists.append(row[0:9].tolist())
line_count = 0

Out_CDS = open(outCDSpath, "w+")

for line in CDS_isoforms_list_of_lists:
    if int(line[3]) == int(line[4]):
        continue
    Out_CDS.writelines('\t'.join(map(str, line)) + "\n")
    line_count += 1
Out_CDS.close()

### In this last part we keep cleaning the file
# Remove genes which are difficult to classify between
# 2 or more genes and genes without a reference.

IsoCDS = pd.read_table(outCDSpath, sep='\t', comment="#", dtype='str', header=None)
IsoCDS.columns = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
Gene_name = IsoCDS["attributes"].str.replace('gene_id "', '', regex=True)
Gene_name = Gene_name.str.replace('"; transcript_id.*', '', regex=True)
IsoCDS["Gene"] = Gene_name
Transcript_name = IsoCDS["attributes"].str.replace('.*transcript_id "', '', regex=True)
Transcript_name = Transcript_name.str.replace('"; .*', '', regex=True)
IsoCDS["Transcript"] = Transcript_name
Gene_name_u = Gene_name.unique().tolist()
IsoCDS["Ref_gene"] = "-"
for index, row in IsoCDS.iterrows():
    ref_ID = len(re.findall("ref_gene_id",row["attributes"]))
    if ref_ID > 0:
        ref_name = re.sub('.*ref_gene_id "', '', row["attributes"])
        ref_name = re.sub('"; .*', '', ref_name)
        row["Ref_gene"] = ref_name

STGTF_onlyIso_transcript = STGTF_onlyIso[:][STGTF_onlyIso.type.isin(["transcript"])].copy()
TranscriptsSTGTF = STGTF_onlyIso_transcript["Transcript"]
Gene_name_iso_u = STGTF_onlyIso_transcript["Gene"].unique().tolist()

TPMs = STGTF_onlyIso_transcript["attributes"].str.replace('.*; TPM "', '', regex=True)
TPMs = TPMs.str.replace('";', '')
STGTF_onlyIso_transcript["TPM"] = TPMs

TPMs_per_transcript = {}
Ref_id_per_transcript = {}

# Create a dictionary to know the TPMs per transcript
# Create a dictionary to indicate if has a reference sequence.
for i in range(0,len(TranscriptsSTGTF)):
    TPMs_per_transcript[TranscriptsSTGTF.tolist()[i]] = float(TPMs.tolist()[i])
    Ref_gene = STGTF_onlyIso_transcript["Ref_gene"].tolist()[i]
    if Ref_gene != "-":
        Ref_id_per_transcript[TranscriptsSTGTF.tolist()[i]] = True
    else:
        Ref_id_per_transcript[TranscriptsSTGTF.tolist()[i]] = False
transcripts_to_remove = []
ncGenes_to_remove = []

# Remove ncRNA
STncGenes = IsoCDS["Gene"][IsoCDS.Ref_gene.isin(ncRNA_Genes.to_list())].to_list()
ncGenes_to_remove.extend(STncGenes)

for gene in Gene_name_iso_u:
    # Detect in one specific gene the subset of transcripts
    IsoCDS_subset = IsoCDS[:][IsoCDS.Gene.isin([gene])]

    transcripts_dict = {}
    Unique_transcripts = IsoCDS_subset["Transcript"].unique().tolist()
    # In this loop we will add the CDS starts and ends creating a unique code for each transcript in the gene
    # If 2 or more transcripts have the same code this means that they have the same CDSs and we only will keep one of them
    for transcript in Unique_transcripts:
        Isoform_CDS_table_sub_subset = IsoCDS_subset[:][IsoCDS_subset.Transcript.isin([transcript])]
        CDS_starts = Isoform_CDS_table_sub_subset["start"].tolist()
        CDS_ends = Isoform_CDS_table_sub_subset["end"].tolist()
        transcripts_dict[transcript] = str("".join(CDS_starts))+str("".join(CDS_ends))

    # Now check if 2 or more transcripts has the same CDS sequence
    # To do se we will create a reverse dictionary and check if there is more than one CDS code repeated

    rev_multidict = {}
    for key, value in transcripts_dict.items():
        rev_multidict.setdefault(value, set()).add(key)
    Repeated_transcripts_set_list = [values for key, values in rev_multidict.items() if len(values) > 1]

    if len(Repeated_transcripts_set_list) > 0:
        Repeated_transcripts = list(Repeated_transcripts_set_list[0])

        # Now check it some of the transcripts has a ref ID
        Repeated_transcripts_refID_subset = {}
        Repeated_transcripts_TPMs_subset = {}
        TF_refID_list = []
        TPM_list = []
        for i in Repeated_transcripts:
            Repeated_transcripts_refID_subset[i] = Ref_id_per_transcript[i]
            TF_refID_list.append(Ref_id_per_transcript[i])
            Repeated_transcripts_TPMs_subset[i] = TPMs_per_transcript[i]
            TPM_list.append(TPMs_per_transcript[i])
        # Now we have 2 outcomes: 1 One of the sequences have a reference ID, then remove all except that one
        # 2 None of them have a ref ID so we have to check the best option with TPMs
        TPM_list.sort()
        if True in TF_refID_list:
            list_of_transcripts_to_remove = {key:value for key, val in Repeated_transcripts_refID_subset.items() if val != True}
            transcripts_to_remove.extend(list(list_of_transcripts_to_remove))
        else:
            max_TPM = TPM_list[-1]
            list_of_transcripts_to_remove = {key:value for key, val in Repeated_transcripts_TPMs_subset.items() if val < max_TPM}
            transcripts_to_remove.extend(list(list_of_transcripts_to_remove))

All_transcripts_with_CDS = IsoCDS["Transcript"].unique().tolist()
Transcripts_with_CDS_kept = list(set(All_transcripts_with_CDS).difference(transcripts_to_remove))
All_genes_with_CDS = IsoCDS["Gene"].unique().tolist()
Genes_with_CDS_kept = list(set(All_genes_with_CDS).difference(ncGenes_to_remove))

Isoforms_transcript_table_tsv = STGTF_onlyIso_transcript[['Gene', 'Transcript', 'Ref_gene', "TPM"]]
Isoforms_transcript_table_tsv_kept = Isoforms_transcript_table_tsv[:][Isoforms_transcript_table_tsv.Transcript.isin(Transcripts_with_CDS_kept)]
Isoforms_transcript_table_tsv_kept_v2 = Isoforms_transcript_table_tsv_kept[:][Isoforms_transcript_table_tsv_kept.Gene.isin(Genes_with_CDS_kept)]

Genes = Isoforms_transcript_table_tsv_kept_v2["Gene"].tolist()
Transcripts_per_gene = {i: Genes.count(i) for i in Genes}
Genes_with_isoforms = []
for key, value in Transcripts_per_gene.items():
    if value > 1:
        Genes_with_isoforms.append(key)

Isoforms_transcript_table_tsv_kept_still_isoforms = Isoforms_transcript_table_tsv_kept_v2[:][Isoforms_transcript_table_tsv_kept_v2.Gene.isin(Genes_with_isoforms)]
# pd.DataFrame(Isoforms_transcript_table_tsv_kept_still_isoforms).to_csv(Isoform_tsv_path, sep = "\t", index = False, header = True, quoting= csv.QUOTE_NONE)

Good_isoforms = STGTF_onlyIso[:][STGTF_onlyIso.Transcript.isin(Isoforms_transcript_table_tsv_kept_still_isoforms["Transcript"].tolist())]
Good_isoforms = Good_isoforms.drop(columns = "Transcript")

pd.DataFrame(Good_isoforms).to_csv(IsoTranscriptsPath, sep = "\t", index = False, header = False, quoting= csv.QUOTE_NONE)

Good_isoforms_CDS = IsoCDS[:][IsoCDS.Transcript.isin(Isoforms_transcript_table_tsv_kept_still_isoforms["Transcript"].tolist())]
Good_isoforms_CDS = Good_isoforms_CDS.drop(columns = "Transcript")
Good_isoforms_CDS = Good_isoforms_CDS.drop(columns = "Gene")
Good_isoforms_CDS = Good_isoforms_CDS.drop(columns = "Ref_gene")

pd.DataFrame(Good_isoforms_CDS).to_csv(IsoCDSsPath, sep = "\t", index = False, header = False, quoting= csv.QUOTE_NONE)
