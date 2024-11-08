import sys

def parse_arguments():
    try:
        _, gff_compare_tracking_path, symbols_chosen, gxf_to_extract_path, gxf_output_path = sys.argv
    except ValueError:
        print("error: gffCompareExtraction.py <gffCompare.tracking file> <gffCompare identifiers chosen separated by ','> <gff to extract sequences> <Output gff>")
        sys.exit(1)

    symbols_chosen_list = symbols_chosen.split(",") if "," in symbols_chosen else [symbols_chosen]
    return gff_compare_tracking_path, symbols_chosen_list, gxf_to_extract_path, gxf_output_path

def extract_desired_symbols(gff_compare_tracking_path, symbols_chosen_list):
    desired_genes, desired_transcripts = [], []
    with open(gff_compare_tracking_path, "r") as file:
        for line in file:
            fields = line.split("\t")
            if fields[3] in symbols_chosen_list:
                desired_genes.append(fields[4].split(":")[1].split("|")[0])
                desired_transcripts.append(fields[4].split(":")[1].split("|")[1])
    return desired_genes, desired_transcripts

def map_transcripts_to_genes(gxf_to_extract_path, desired_transcripts):
    gene_transcripts, multi_transcripts_valids = {}, []
    with open(gxf_to_extract_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            feature, attr = fields[2], fields[8]
            if feature in ["transcript"]:
                transcript_name = attr.split("ID=")[1].split(";")[0]
                gene_name = attr.split("Parent=")[1].split(";")[0]
                if transcript_name in desired_transcripts:
                    gene_transcripts.setdefault(gene_name, []).append(transcript_name)
                    if len(gene_transcripts[gene_name]) > 1:
                        multi_transcripts_valids.append(gene_name)
    return gene_transcripts, multi_transcripts_valids

def count_exons_per_transcript(gxf_to_extract_path, desired_transcripts):
    exons_in_transcript = {}
    with open(gxf_to_extract_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            feature, attr = fields[2], fields[8]
            if feature == "exon":
                transcript_name = attr.split("Parent=")[1].split(";")[0]
                if transcript_name in desired_transcripts:
                    exons_in_transcript[transcript_name] = exons_in_transcript.get(transcript_name, 0) + 1
    return exons_in_transcript

def select_valid_transcripts(gene_transcripts, multi_transcripts_valids, exons_in_transcript):
    valid_transcripts = []
    for gene, transcripts in gene_transcripts.items():
        if gene in multi_transcripts_valids:
            min_exons = float('inf')
            for transcript in transcripts:
                exon_count = exons_in_transcript.get(transcript, float('inf'))
                if exon_count < min_exons:
                    min_exons = exon_count
                    selected_transcript = transcript
            valid_transcripts.append(selected_transcript)
        else:
            valid_transcripts.append(transcripts[0])
    return valid_transcripts

def write_output_gff(gxf_to_extract_path, desired_genes, desired_transcripts, valid_transcripts, gxf_output_path):
    with open(gxf_to_extract_path, "r") as infile, open(gxf_output_path, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            fields = line.split("\t")
            feature, attr = fields[2], fields[8]
            if feature == "gene":
                gene_name = attr.split("ID=")[1].split(";")[0]
                if gene_name in desired_genes:
                    outfile.write(line)
            elif feature in ["transcript"]:
                transcript_name = attr.split("ID=")[1].split(";")[0]
                if transcript_name in desired_transcripts and transcript_name in valid_transcripts:
                    outfile.write(line)
            else:
                transcript_name = attr.split("Parent=")[1].split(";")[0]
                if transcript_name in desired_transcripts and transcript_name in valid_transcripts:
                    outfile.write(line)

def main():
    gff_compare_tracking_path, symbols_chosen_list, gxf_to_extract_path, gxf_output_path = parse_arguments()
    desired_genes, desired_transcripts = extract_desired_symbols(gff_compare_tracking_path, symbols_chosen_list)
    gene_transcripts, multi_transcripts_valids = map_transcripts_to_genes(gxf_to_extract_path, desired_transcripts)
    exons_in_transcript = count_exons_per_transcript(gxf_to_extract_path, desired_transcripts)
    valid_transcripts = select_valid_transcripts(gene_transcripts, multi_transcripts_valids, exons_in_transcript)
    write_output_gff(gxf_to_extract_path, desired_genes, desired_transcripts, valid_transcripts, gxf_output_path)

if __name__ == "__main__":
    main()
