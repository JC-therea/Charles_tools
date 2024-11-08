import sys

def parse_arguments():
    try:
        _, ref_path, original_fun_path, filtered_fun_path, rescued_out_path = sys.argv
    except ValueError:
        print("Usage: gffCompareExtraction.py <reference.gff> <funannotate.original.gff> <funannotate.filtered.gff> <Output gff>")
        sys.exit(1)
    return ref_path, original_fun_path, filtered_fun_path, rescued_out_path

def parse_filtered_file(filtered_fun_path):
    old_genes, old_transcripts = set(), set()
    with open(filtered_fun_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            feature, attr = fields[2], fields[8]
            if feature == "gene":
                gene_name = extract_attribute(attr, "ID")
                old_genes.add(gene_name)
            elif feature in ["mRNA", "rRNA", "ncRNA", "tRNA"]:
                transcript_name = extract_attribute(attr, "ID")
                old_transcripts.add(transcript_name)
            else:
                transcript_name = extract_attribute(attr, "Parent")
                old_transcripts.add(transcript_name)
    return old_genes, old_transcripts

def extract_attribute(attr, key):
    try:
        return next(item.split("=")[1] for item in attr.split(";") if item.startswith(key))
    except StopIteration:
        return None

def parse_transcripts(original_fun_path):
    transcript_to_gene = {}
    with open(original_fun_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] == "mRNA":
                transcript_name = extract_attribute(fields[8], "ID")
                gene_name = extract_attribute(fields[8], "Parent")
                transcript_to_gene[transcript_name] = gene_name
    return transcript_to_gene

def parse_utrs(original_fun_path, transcript_to_gene):
    utr3_dict, utr5_dict = {}, {}
    with open(original_fun_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            feature, start, end = fields[2], int(fields[3]), int(fields[4])
            transcript = extract_attribute(fields[8], "Parent")
            gene = transcript_to_gene.get(transcript)

            if feature == "five_prime_UTR":
                utr5_dict[gene] = update_utr(utr5_dict, gene, start, end)
            elif feature == "three_prime_UTR":
                utr3_dict[gene] = update_utr(utr3_dict, gene, start, end)
    return utr3_dict, utr5_dict

def update_utr(utr_dict, gene, start, end):
    if gene not in utr_dict:
        utr_dict[gene] = [start, end]
    else:
        utr_dict[gene][0] = min(utr_dict[gene][0], start)
        utr_dict[gene][1] = max(utr_dict[gene][1], end)
    return utr_dict[gene]

def parse_reference(ref_path, old_genes, transcript_to_gene_ref):
    exon_storage, transcripts_with_cds, genes_with_cds = {}, set(), set()
    with open(ref_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            feature, attr = fields[2], fields[8]
            ref_transcript = extract_attribute(attr, "Parent")
            if feature == "mRNA":
                ref_gene = extract_attribute(attr, "Parent")
                transcript_to_gene_ref[ref_transcript] = ref_gene
            elif feature == "exon":
                ref_gene = transcript_to_gene_ref.get(ref_transcript)
                if ref_gene and ref_gene not in old_genes:
                    exon_storage.setdefault(ref_gene, []).extend([fields[3], fields[4]])
            elif feature == "CDS":
                ref_gene = transcript_to_gene_ref.get(ref_transcript)
                transcripts_with_cds.add(ref_transcript)
                if ref_gene:
                    genes_with_cds.add(ref_gene)
    return exon_storage, transcripts_with_cds, genes_with_cds

def rescue_features(ref_path, old_genes, utr3_dict, utr5_dict, genes_with_cds, transcripts_with_cds, exon_storage, transcript_to_gene_ref, gxf_output):
    exon_counting = {}
    with open(ref_path, "r") as file, open(gxf_output, "w") as out_file:
        for line in file:
            if line.startswith("#"):
                out_file.write(line)
                continue
            fields = line.split("\t")
            feature, start, end, strand, attr = fields[2], int(fields[3]), int(fields[4]), fields[6], fields[8]

            if feature == "gene":
                ref_gene = extract_attribute(attr, "ID")
                if ref_gene not in old_genes and ref_gene in genes_with_cds:
                    start, end = update_positions_with_utrs(ref_gene, strand, start, end, utr5_dict, utr3_dict)
                    out_file.write(format_gff_line(fields, ref_gene, start, end))
            elif feature == "mRNA":
                ref_gene = extract_attribute(attr, "Parent")
                ref_transcript = extract_attribute(attr, "ID")
                if ref_gene not in old_genes and ref_transcript in transcripts_with_cds:
                    start, end = update_positions_with_utrs(ref_gene, strand, start, end, utr5_dict, utr3_dict)
                    out_file.write(format_gff_line(fields, ref_gene + "-T1", start, end, ref_gene))
            elif feature == "exon":
                ref_transcript = extract_attribute(attr, "Parent")
                ref_gene = transcript_to_gene_ref.get(ref_transcript)
                if ref_gene not in old_genes and ref_transcript in transcripts_with_cds:
                    if len(exon_storage[ref_gene]) > 2:
                        start, end = update_positions_for_multiexon(ref_gene, strand, start, end, exon_storage, utr5_dict, utr3_dict)
                    exon_counting[ref_gene] = exon_counting.get(ref_gene, 0) + 1
                    out_file.write(format_gff_line(fields, ref_gene + "-T1.exon" + str(exon_counting[ref_gene]), start, end, ref_gene + "-T1"))
            elif feature == "CDS":
                ref_transcript = extract_attribute(attr, "Parent")
                ref_gene = transcript_to_gene_ref.get(ref_transcript)
                if ref_gene:
                    out_file.write(format_gff_line(fields, ref_gene + "-T1.cds", start, end, ref_gene + "-T1"))

def update_positions_with_utrs(ref_gene, strand, start, end, utr5_dict, utr3_dict):
    if strand == "+":
        start = min(start, int(utr5_dict.get(ref_gene, [start])[0]))
        end = max(end, int(utr3_dict.get(ref_gene, [end])[1]))
    else:
        start = min(start, int(utr3_dict.get(ref_gene, [start])[0]))
        end = max(end, int(utr5_dict.get(ref_gene, [end])[1]))
    return start, end

def update_positions_for_multiexon(ref_gene, strand, start, end, exon_storage, utr5_dict, utr3_dict):
    min_exon = min(exon_storage[ref_gene])
    max_exon = max(exon_storage[ref_gene])
    if strand == "+":
        if start == min_exon:
            start = min(start, int(utr5_dict.get(ref_gene, [start])[0]))
        if end == max_exon:
            end = max(end, int(utr3_dict.get(ref_gene, [end])[1]))
    else:
        if start == min_exon:
            start = min(start, int(utr3_dict.get(ref_gene, [start])[0]))
        if end == max_exon:
            end = max(end, int(utr5_dict.get(ref_gene, [end])[1]))
    return start, end

def format_gff_line(fields, new_id, start, end, parent=None):
    attr = f"ID={new_id};" + (f"Parent={parent};" if parent else "") + "\n"
    fields[3], fields[4], fields[8] = str(start), str(end), attr
    return "\t".join(fields)

def main():
    ref_path, original_fun_path, filtered_fun_path, rescued_out_path = parse_arguments()
    old_genes, old_transcripts = parse_filtered_file(filtered_fun_path)
    transcript_to_gene = parse_transcripts(original_fun_path)
    utr3_dict, utr5_dict = parse_utrs(original_fun_path, transcript_to_gene)
    transcript_to_gene_ref = {}
    exon_storage, transcripts_with_cds, genes_with_cds = parse_reference(ref_path, old_genes, transcript_to_gene_ref)
    rescue_features(ref_path, old_genes, utr3_dict, utr5_dict, genes_with_cds, transcripts_with_cds, exon_storage, transcript_to_gene_ref, rescued_out_path)

if __name__ == "__main__":
    main()
