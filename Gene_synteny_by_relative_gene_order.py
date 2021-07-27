import re
import math
import sys
from datetime import datetime

setting_file_path = sys.argv[1]

try:
    setting_file = open(setting_file_path, "r").readlines()

    for line in setting_file:
        if line.startswith("Reference fasta sequences:"):
            path_ref_fasta = line.replace("\n", "").replace("Reference fasta sequences:","")
            continue
        if line.startswith("Outgroup fasta sequences:"):
            path_orth_fasta = line.replace("\n", "").replace("Outgroup fasta sequences:","")
            continue
        if line.startswith("Reference gff annotation file:"):
            path_ref_gff = line.replace("\n", "").replace("Reference gff annotation file:","")
            continue
        if line.startswith("Outgroup gff annotation file:"):
            path_orth_gff = line.replace("\n", "").replace("Outgroup gff annotation file:","")
            continue
        if line.startswith("Reference and orthologue orthofinder output:"):
            path_orthologues = line.replace("\n", "").replace("Reference and orthologue orthofinder output:","")
            continue
        if line.startswith("Gene-B phylo age path:"):
            GeneB_age_path = line.replace("\n", "").replace("Gene-B phylo age path:","")
            continue
        if line.startswith("Nodes:"):
            Desired_nodes = line.replace("\n", "").replace("Nodes:","")
            continue
        if line.startswith("Check only fasta genes?:"):
            check_fasta_orth = line.replace("\n", "").replace("Check only fasta genes?:","")
            continue
        if line.startswith("Output path:"):
            path_output = line.replace("\n", "").replace("Output path:","")
            continue




except:
    setting_file = open(setting_file_path, "w+")
    setting_file.writelines("Reference fasta sequences:" + "\n")
    setting_file.writelines("Outgroup fasta sequences:" + "\n")
    setting_file.writelines("Reference gff annotation file:" + "\n")
    setting_file.writelines("Outgroup gff annotation file:" + "\n")
    setting_file.writelines("Reference and orthologue orthofinder output:" + "\n")
    setting_file.writelines("Gene-B phylo age path:" + "\n")
    setting_file.writelines("Nodes:" + "\n")
    setting_file.writelines("Check only fasta genes?:" + "\n")
    setting_file.writelines("Output path:" + "\n")
    print("\n" + "Settings file introduced do not exist, creating new settings file with users input...")
    print("\n" + "File created please fill the information")
    exit()

try:
    open(path_ref_fasta,"r")
except:
    print("\n" + "Reference fasta sequences not found")
    exit()
try:
    open(path_ref_gff,"r")
except:
    print("\n" + "Reference gff annotation file not found")
    exit()
try:
    open(path_orth_fasta,"r")
except:
    print("\n" + "Outgroup fasta sequences not found")
    exit()
try:
    open(path_orth_gff,"r")
except:
    print("\n" + "Outgroup gff annotation file not found")
    exit()
try:
    open(path_orthologues,"r")
except:
    print("\n" + "Orthologues file not found")
    exit()
try:
    open(GeneB_age_path,"r")
except:
    print("\n" + "Gene-B phylo file not found")
    exit()


# Function to get the names that we want to include

def get_header(filename):
    with open(filename, 'r') as fd:
        headers = []
        for line in fd.readlines():
            if '>' in line:
                headers.append(line.strip()[1:])
            else:
                pass
    return headers

# Function to get positions

def genes_pos(filename, admitted_genes, only_fasta_IDs):
    with open(filename, 'r') as fd:
        position = 1
        genes_to_pos = {}
        pos_to_genes = {}
        terminal_genelist = []
        gene_to_chr = {}
        current_chr=""
        for line in fd.readlines():
            line_ref_without_n = line.replace("\n", "")
            Chr, Source, Feature, Start, End, Score, Strand, Frame, Attributes = re.split("\t", line_ref_without_n)

            if Feature == "mRNA":
                    Gene_ID = re.sub(r"ID=","",re.sub(r';.*','',Attributes))

                    if only_fasta_IDs == "True":
                        if Gene_ID in admitted_genes:
                            gene_to_chr[Gene_ID] = Chr
                            genes_to_pos[Gene_ID] = position
                            pos_to_genes[position] = Gene_ID
                            position += 1
                        else:
                            continue
                    else:
                        gene_to_chr[Gene_ID] = Chr
                        genes_to_pos[Gene_ID] = position
                        pos_to_genes[position] = Gene_ID
                        position += 1

                    if current_chr != Chr:
                        if current_chr == "":
                            terminal_genelist.append(Gene_ID)
                            current_chr = Chr
                        else:
                            terminal_genelist.append(Gene_ID)
                            if position - 2 > 0:
                                terminal_genelist.append(pos_to_genes[position - 2])
                                current_chr = Chr

    return genes_to_pos, pos_to_genes, terminal_genelist[:-1], gene_to_chr

list_of_genes_ref = get_header(path_ref_fasta)
list_of_genes_orth = get_header(path_orth_fasta)

genes_to_pos_ref, pos_to_genes_ref, terminal_genelist_ref ,gene_to_chr_ref = genes_pos(path_ref_gff,list_of_genes_ref, "False")
genes_to_pos_orth, pos_to_genes_orth, terminal_genelist_orth ,gene_to_chr_orth = genes_pos(path_orth_gff,list_of_genes_orth, check_fasta_orth)

# Get Orthologues

with open(path_orthologues, 'r') as fd:
    ref_to_orth = {}
    for line in fd.readlines():
        line_ref_without_n = line.replace("\n", "")
        Orthogroup, ref_genes, orthologues = re.split("\t", line_ref_without_n)
        if Orthogroup == "Orthogroup":
            pass
        else:
            ref_genes_list = re.split(", ",ref_genes)
            for gene in ref_genes_list:
                ref_to_orth[gene] = re.split(", ",orthologues)

output_path = path_output + "/Gene_study_position_" + datetime.today().strftime('%Y-%m-%d') + ".tsv"
output = open(output_path,"w+")
first_line = "Gene" +"\t"+ "Closest_left_gene_w_Orth" +"\t"+ "Closest_right_gene_w_Orth" +"\t"+ \
             "Correct_left_orth" +"\t"+ "Correct_right_orth" +"\t"+ "Diff_ref_neighbours" +"\t"+ "Diff_orth" +"\t"+ "Classification" +"\n"
output.write(first_line)

# Read Gene-B Phylo file

with open(GeneB_age_path, 'r') as fd:
    Gene_2_Node = {}
    for line in fd.readlines():
        if line.startswith("Gene"):
            continue
        line_ref_without_n = line.replace("\n", "")
        Gene, Orthogroup, Clade, Node = re.split("\t", line_ref_without_n)
        Gene_2_Node[Gene] = Node

# Until here I have: orthologues "ref_to_orth",
# gene relative position in the reference species "genes_to_pos_ref", "pos_to_genes_ref",
# and relative position in the orthologue species "genes_to_pos_orth", "pos_to_genes_orth"

for position in range(1, len(list_of_genes_ref)):
    studied_gene = pos_to_genes_ref[position]
    left_ref_gene = ""
    right_ref_gene = ""
    left_pos_ref = 0
    right_pos_ref = 0
    if studied_gene not in list_of_genes_ref:
        continue
    if studied_gene in Gene_2_Node.keys():
        if Gene_2_Node[studied_gene] not in Desired_nodes:
            continue

    if studied_gene in terminal_genelist_ref:
        line = studied_gene +"\t"+ "None" +"\t"+ "None" +"\t"+ "None" +"\t"+ "None" +"\t"+ "None" +"\t"+ "None" +"\t"+ "terminal" +"\n"
        output.write(line)
        continue

    while left_ref_gene == "":
        left_pos_ref += (-1)
        left_neighbour_candidate = pos_to_genes_ref[position + left_pos_ref]
        if left_neighbour_candidate in ref_to_orth.keys():
            orth_left = ref_to_orth[left_neighbour_candidate]
            left_ref_gene = left_neighbour_candidate
        else:
            if left_neighbour_candidate in terminal_genelist_ref:
                left_ref_gene = "terminal no orthologue"
                break
            else:
                continue

    while right_ref_gene == "":
        right_pos_ref += (+1)
        right_neighbour_candidate = pos_to_genes_ref[position + right_pos_ref]
        if right_neighbour_candidate in ref_to_orth.keys():
            orth_right = ref_to_orth[right_neighbour_candidate]
            right_ref_gene = right_neighbour_candidate
        else:
            if right_neighbour_candidate in terminal_genelist_ref:
                right_ref_gene = "terminal no orthologue"
                break
            else:
                continue

    if left_ref_gene == "terminal no orthologue" or right_ref_gene == "terminal no orthologue":
        line = studied_gene + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "terminal neighbours" + "\n"
        output.write(line)
        continue

    dist_ref = abs(right_pos_ref - left_pos_ref)
    list_orth_distances = []

    for R in orth_right:
        list_orth_distances.append([abs(genes_to_pos_orth[R] - genes_to_pos_orth[L]) for L in orth_left])

    final_distance_list_orth = [item for t in list_orth_distances for item in t]
    orthologues_distance = min(final_distance_list_orth)
    min_index = final_distance_list_orth.index(orthologues_distance)
    orth_right_OK = orth_right[math.floor(min_index/len(orth_left))]

    if (min_index % len(orth_left)) == 0:
        orth_left_OK = orth_left[0]
    else:
        orth_left_OK = orth_left[(min_index % len(orth_left)) - 1]

    reference_distance = genes_to_pos_ref[right_ref_gene] - genes_to_pos_ref[left_ref_gene]
    if orthologues_distance == 0:
        line = studied_gene + "\t" + left_ref_gene + "\t" + right_ref_gene + "\t" + orth_left_OK + "\t" + orth_right_OK + "\t" + str(
            reference_distance) + "\t" + str(orthologues_distance) + "\t" + "Duplicated neighbours" + "\n"
        output.write(line)
    elif gene_to_chr_orth[orth_left_OK] != gene_to_chr_orth[orth_right_OK]:
        line = studied_gene + "\t" + left_ref_gene + "\t" + right_ref_gene + "\t" + orth_left_OK + "\t" + orth_right_OK + "\t" + str(
            reference_distance) + "\t" + str(orthologues_distance) + "\t" + "Genomic rearrengement" + "\n"
        output.write(line)
    elif orthologues_distance < reference_distance:
        if abs(orthologues_distance) == 1:
            line = studied_gene + "\t" + left_ref_gene + "\t" + right_ref_gene + "\t" + orth_left_OK + "\t" + orth_right_OK + "\t" + str(reference_distance) + "\t" + str(orthologues_distance) + "\t" + "De novo" + "\n"
            output.write(line)
        else:
            line = studied_gene + "\t" + left_ref_gene + "\t" + right_ref_gene + "\t" + orth_left_OK + "\t" + orth_right_OK + "\t" + str(reference_distance) + "\t" + str(orthologues_distance) + "\t" + "Possible_de_novo" + "\n"
            output.write(line)
    elif orthologues_distance == reference_distance:
        line = studied_gene + "\t" + left_ref_gene + "\t" + right_ref_gene + "\t" + orth_left_OK + "\t" + orth_right_OK + "\t" + str(
            reference_distance) + "\t" + str(orthologues_distance) + "\t" + "Divergence" + "\n"
        output.write(line)
    else:
        line = studied_gene + "\t" + left_ref_gene + "\t" + right_ref_gene + "\t" + orth_left_OK + "\t" + orth_right_OK + "\t" + str(
            reference_distance) + "\t" + str(orthologues_distance) + "\t" + "Gene inclusion" + "\n"
        output.write(line)

output.close()

### Print the position of each species
out_positions_ref_path = path_output + "/Reference_gene_position_" + datetime.today().strftime('%Y-%m-%d') + ".tsv"
out_positions_ref = open(out_positions_ref_path,"w+")
for position in range(1, len(list_of_genes_ref)):
    studied_gene = pos_to_genes_ref[position]
    out_positions_ref.writelines(studied_gene + "\t" + str(position) + "\n")
out_positions_ref.close()

out_positions_outgroup_path = path_output + "/Outgroup_gene_position_" + datetime.today().strftime('%Y-%m-%d') + ".tsv"
out_positions_outgroup = open(out_positions_outgroup_path,"w+")
for position in range(1, len(pos_to_genes_orth.keys())):
    studied_gene = pos_to_genes_orth[position]
    out_positions_outgroup.writelines(studied_gene + "\t" + str(position) + "\n")

out_positions_outgroup.close()
