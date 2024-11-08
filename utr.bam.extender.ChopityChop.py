import sys

def parse_arguments():
    try:
        _, ref_path, rescued_out_path, chopped_utrs_path = sys.argv
    except ValueError:
        print("Usage: ChopityChop.py <reference.gff> <funannotate.rescued.gff> <Output gff>")
        sys.exit(1)
    return ref_path, rescued_out_path, chopped_utrs_path

def load_features(input_file, feature_types):
    """Parse GFF file and return a dictionary with specified feature types and their coordinates."""
    features_dict = {}
    with open(input_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chr, strand, feature, start, end = fields[0], fields[6], fields[2], int(fields[3]), int(fields[4])
            attr = fields[8]
            if feature in feature_types:
                features_dict.setdefault((chr, strand), []).append((start, end))
    return features_dict

def parse_transcripts(rescued_out_path):
    transcripts_to_genes, genes_to_transcripts, transcripts_with_utrs = {}, {}, []
    with open(rescued_out_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            feature, attr = fields[2], fields[8]
            if feature in ["three_prime_UTR", "five_prime_UTR"]:
                transcripts_with_utrs.append(get_attribute(attr, "Parent"))
            if feature in ["transcript"]:
                transcript_name = get_attribute(attr, "ID")
                gene_name = get_attribute(attr, "Parent")
                transcripts_to_genes[transcript_name] = gene_name
                genes_to_transcripts[gene_name] = transcript_name
    return transcripts_to_genes, genes_to_transcripts, set(transcripts_with_utrs)

def get_attribute(attr, key):
    try:
        return next(item.split("=")[1] for item in attr.split(";") if item.startswith(key))
    except StopIteration:
        return None

def adjust_utr_positions(line, cds_trna_dict, transcripts_with_utrs):
    fields = line.strip().split("\t")
    chr, strand, feature, start, end = fields[0], fields[6], fields[2], int(fields[3]), int(fields[4])
    attr = fields[8]
    transcript_name = get_attribute(attr, "Parent" if feature in ["exon", "three_prime_UTR", "five_prime_UTR"] else "ID")

    if transcript_name not in transcripts_with_utrs:
        return line

    for cds_start, cds_end in cds_trna_dict.get((chr, strand), []):
        if start == cds_start and end == cds_end or attr.startswith("ID=gene-YNC"):
            continue
        if start < cds_end and end > cds_start:
            left_rest, right_rest = cds_start - start, end - cds_end
            if strand == "+":
                if feature == "five_prime_UTR":
                    start = cds_end + 1
                elif feature == "three_prime_UTR":
                    end = cds_start - 1
                elif left_rest < right_rest:
                    start = cds_end + 1
                else:
                    end = cds_start - 1
            else:
                if feature == "five_prime_UTR":
                    start = cds_start - 1
                elif feature == "three_prime_UTR":
                    end = cds_end + 1
                elif left_rest < right_rest:
                    start = cds_end + 1
                else:
                    end = cds_start - 1
    return "\t".join(fields[:3] + [str(start), str(end)] + fields[5:])

def chop_utrs(reference_path, rescued_out_path, chopped_utrs_path, cds_trna_dict, transcripts_with_utrs):
    with open(reference_path, "r") as input_file, open(chopped_utrs_path, "w") as output_file:
        for line in input_file:
            if line.startswith("#"):
                output_file.write(line)
                continue
            fields = line.strip().split("\t")
            feature = fields[2]
            if feature in ["three_prime_UTR", "five_prime_UTR", "gene", "transcript", "exon"]:
                modified_line = adjust_utr_positions(line, cds_trna_dict, transcripts_with_utrs)
                output_file.write(modified_line + "\n")
            else:
                output_file.write(line + "\n")

def main():
    ref_path, rescued_out_path, chopped_utrs_path = parse_arguments()
    cds_trna_dict = load_features(ref_path, ["tRNA"])
    transcripts_to_genes, genes_to_transcripts, transcripts_with_utrs = parse_transcripts(rescued_out_path)
    chop_utrs(ref_path, rescued_out_path, chopped_utrs_path, cds_trna_dict, transcripts_with_utrs)

if __name__ == "__main__":
    main()
