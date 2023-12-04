import re
from Bio.Seq import Seq
from Bio import SeqIO
import math
import pandas as pd
from datetime import datetime


##########################################################################
#PATHS
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Report the number of times each triplet appears in a given multi fasta file')
    parser.add_argument('-g', '--genome', required=True, help='Path to the genome in fasta format')
    parser.add_argument('-gff', '--annotation_file', required=True, help='Path to the annotation file in gff format')
    parser.add_argument('-vcf', '--vcf', required=True, help='Path to the annotation file in vcf format')
    parser.add_argument('-o', '--output_file', required=True, help='Output file in tsv format')
    return parser.parse_args()

args = get_args()
outputPath = args.output_file + "PNPS_per_protein.tsv"
#CDSPath = args.output_file + "CDS_N_new_method.fa"

GENOME = args.genome
GFF = args.annotation_file
VCF_path = args.vcf
Output_path = outputPath

Output = open(outputPath,"w+")
#written_sequences = open(CDSPath, "w+")


def create_codons(Nucleotide_sequence):
    seq_l = len(Nucleotide_sequence)
    codons = []
    start = 0
    while start < seq_l:
        codons.append([Nucleotide_sequence[start:start + 3]])
        start += 3
    return (codons)


Output.write("Protein" + "\t" + "PN" + "\t" + "PS" + "\t" + "Fail_references" + "\n")


def PNPS_searcher(sequence, start_exons, end_exons, SNP_positions, SNP_alternative_nucleotide, SNP_reference_nucleotide,
                  Strand):
    SNPs_local_pos = []
    SNPs_local_alt = []
    SNPs_local_ref = []
    full_CDS = []
    accumulated_length = 0
    PN = 0
    PS = 0

    Fail_reference = 0

    if Strand == "+":

        for i in range(0, len(start_exons)):

            global_exon_start = start_exons[i]
            global_exon_end = end_exons[i]
            exon = sequence[global_exon_start:global_exon_end]

            SNPs_pos_in_exon = SNP_positions[(SNP_positions >= global_exon_start) & (SNP_positions < global_exon_end)]
            SNPs_alt_in_exon = SNP_alternative_nucleotide[
                (SNP_positions >= global_exon_start) & (SNP_positions < global_exon_end)]
            SNPs_ref_in_exon = SNP_reference_nucleotide[
                (SNP_positions >= global_exon_start) & (SNP_positions < global_exon_end)]

            full_CDS.append(exon)
            SNP_local_positions = SNPs_pos_in_exon - global_exon_start + accumulated_length

            for h in range(0, len(SNP_local_positions)):
                SNPs_local_pos.append(int(pd.array(SNP_local_positions)[h]))
                SNPs_local_alt.append(pd.array(SNPs_alt_in_exon)[h])
                SNPs_local_ref.append(pd.array(SNPs_ref_in_exon)[h])

            full_CDS_seq = ''.join(full_CDS)
            accumulated_length += len(exon)

        Affected_codons_aux = [math.floor(round_number / 3) for round_number in SNPs_local_pos]
        Affected_codons = pd.array(Affected_codons_aux)

        for j in range(0, len(SNPs_local_pos)):

            SNP_pos = SNPs_local_pos[j]
            Alt = SNPs_local_alt[j]
            Ref = SNPs_local_ref[j]
            Alt_seq = full_CDS_seq[:SNP_pos] + Alt + full_CDS_seq[SNP_pos + 1:]
            Affected_codon = Affected_codons[j]
            # print(SNP_pos)
            if Ref != full_CDS_seq[SNP_pos]:
                Fail_reference += 1

            ref_codon = ''.join(create_codons(full_CDS_seq)[Affected_codon])
            alt_codon = ''.join(create_codons(Alt_seq)[Affected_codon])

            AA_ref = str(Seq(ref_codon).translate())
            AA_alt = str(Seq(alt_codon).translate())

            if AA_ref == AA_alt:
                PS += 1
            else:
                PN += 1
        if len(SNPs_local_pos) > 0:
            fail_prop = float(Fail_reference) / len(SNPs_local_pos) * 100
        else:
            fail_prop = 0

        #written_sequences.write(full_CDS_seq + "\n")

        return (PN, PS, fail_prop)

    elif Strand == '-':

        for i in range(len(start_exons) - 1, -1, -1):

            global_exon_start = start_exons[i]
            global_exon_end = end_exons[i]
            exon = sequence[global_exon_start:global_exon_end]

            exon_seq = str(Seq(exon).reverse_complement())
            full_CDS.append(exon_seq)

            SNPs_pos_in_exon = SNP_positions[(SNP_positions >= global_exon_start) & (SNP_positions < global_exon_end)]
            SNPs_alt_in_exon = SNP_alternative_nucleotide[
                (SNP_positions >= global_exon_start) & (SNP_positions < global_exon_end)]
            SNPs_ref_in_exon = SNP_reference_nucleotide[
                (SNP_positions >= global_exon_start) & (SNP_positions < global_exon_end)]

            SNP_local_positions = global_exon_end - SNPs_pos_in_exon - 1 + accumulated_length

            for h in range(0, len(SNP_local_positions)):
                SNPs_local_pos.append(int(pd.array(SNP_local_positions)[h]))
                SNPs_local_alt.append(pd.array(SNPs_alt_in_exon)[h])
                SNPs_local_ref.append(pd.array(SNPs_ref_in_exon)[h])

            full_CDS_seq = ''.join(full_CDS)
            accumulated_length += len(exon)

        Affected_codons_aux = [math.floor(round_number / 3) for round_number in SNPs_local_pos]
        Affected_codons = pd.array(Affected_codons_aux)

        for j in range(0, len(SNPs_local_pos)):

            SNP_pos = SNPs_local_pos[j]
            Alt = str(Seq(SNPs_local_alt[j]).complement())
            Ref = str(Seq(SNPs_local_ref[j]).complement())
            # print(SNP_pos)
            if Ref != full_CDS_seq[SNP_pos]:
                Fail_reference += 1

            Alt_seq = full_CDS_seq[:SNP_pos] + Alt + full_CDS_seq[SNP_pos + 1:]
            Affected_codon = Affected_codons[j]

            ref_codon = ''.join(create_codons(full_CDS_seq)[Affected_codon])
            alt_codon = ''.join(create_codons(Alt_seq)[Affected_codon])

            AA_ref = str(Seq(ref_codon).translate())
            AA_alt = str(Seq(alt_codon).translate())

            if AA_ref == AA_alt:
                PS += 1
            else:
                PN += 1
        if len(SNPs_local_pos) > 0:
            fail_prop = float(Fail_reference) / len(SNPs_local_pos) * 100
        else:
            fail_prop = 0
        #written_sequences.write(full_CDS_seq + "\n")
        return (PN, PS, fail_prop)


############################################################################
############################################################################

# To know when the script starts

globalStartTime = datetime.now()

#######################################################
#######################################################
# Load VCF file in Pandas data frame

print("Loading VCF file...")

startTime = datetime.now()

SNP_table = pd.read_table(VCF_path, sep='\t', comment="#", dtype='str', header=None)
print("Transform column 2 into numeric...")
SNP_table_v2 = SNP_table.rename(columns={0: "Chr", 1: "Position", 3: "Reference", 4: "Alternative"})
SNP_table_v2["Position"] = pd.to_numeric(SNP_table_v2["Position"])

print("Load VCF file in " + str(datetime.now() - startTime))

#######################################################
#######################################################
# Load annotation table in pandas

print("Loading annotation file...")
startTime = datetime.now()

CDS_table = pd.read_table(GFF, sep='\t', comment="#", dtype='str', header=None)
CDS_table_v2 = CDS_table.rename(columns={0: "Chr", 3: "Start", 4: "End", 6: "Strand", 8: "Protein"})
CDS_table_v2.Start = pd.to_numeric(CDS_table_v2.Start)
CDS_table_v2.End = pd.to_numeric(CDS_table_v2.End)
print("Load GFF file in " + str(datetime.now() - startTime))

#######################################################
#######################################################

print("Reading genome by chromosomes...")
startTime = datetime.now()

# Change this part to make it variable from the vcf file
allowed_chr = ["2L", "2R", "3L", "3R", "4", "X"]
Proteins_PNPS = {}

for protein in CDS_table_v2.Protein.unique():
    Proteins_PNPS[protein] = [-1, -1, -1]

fasta_sequences = SeqIO.parse(open(GENOME), 'fasta')

for fasta in fasta_sequences:
    name = fasta.id.split(" ")[0]
    # print(name)
    if name in allowed_chr:
        CDS_table_v2_subset = CDS_table_v2[:][CDS_table_v2.Chr.isin([name])]
        SNP_table_v2_subset = SNP_table_v2[:][SNP_table_v2.Chr.isin([name])]

        chromosome_seq = str(fasta.seq)

        Prots_in_chr = CDS_table_v2_subset.Protein.unique()

        for Protein in Prots_in_chr:
            CDS_table_v2_protein_subset = CDS_table_v2_subset[:][CDS_table_v2_subset.Protein.isin([Protein])]
            Strand = CDS_table_v2_protein_subset.Strand.values[0]
            Starts = CDS_table_v2_protein_subset.Start.values
            Ends = CDS_table_v2_protein_subset.End.values
            # print(Protein)
            # print(Strand)
            #written_sequences.write(">" + Protein + "\n")
            Proteins_PNPS[Protein][0], Proteins_PNPS[Protein][1], Proteins_PNPS[Protein][2] = PNPS_searcher(
                chromosome_seq, (Starts - 1), Ends, (SNP_table_v2_subset.Position - 1), SNP_table_v2_subset.Alternative,
                SNP_table_v2_subset.Reference, Strand)

            Output.write(
                Protein + "\t" + str(Proteins_PNPS[Protein][0]) + "\t" + str(Proteins_PNPS[Protein][1]) + "\t" + str(
                    Proteins_PNPS[Protein][2]) + "\n")
        print("Chromosome " + name + ' done in ' + str(datetime.now() - startTime))
    else:
        pass

Output.close()

print("Total script time: " + str(datetime.now() - globalStartTime))
