import sys

def modify_gff(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        exon_count = {}
        cds_count = {}

        for line in infile:
            if line.startswith("#") or line.strip() == "":
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            seqname, source, feature, start, end, score, strand, frame, attributes = fields

            if source != "merged.genePred":
                continue

            attr_dict = parse_attributes(attributes)

            if feature == "gene":
                new_attr = {k: v for k, v in attr_dict.items() if k in ["ID", "gene_id", "Name", "locus_tag"]}
                new_line = format_gff_line(fields, new_attr)
                outfile.write(new_line)

            elif feature == "transcript":
                new_attr = {k: v for k, v in attr_dict.items() if k in ["ID", "Parent", "gene_id"]}

                if new_attr["ID"] == new_attr["Parent"]:
                    new_attr["ID"] = new_attr["ID"] + "-T1"
                    new_attr["transcript_id"] = new_attr["ID"]

                new_line = format_gff_line(fields, new_attr, feature)
                outfile.write(new_line)

            elif feature == "exon":
                new_attr = {k: v for k, v in attr_dict.items() if k in ["ID", "Parent", "gene_id", "transcript_id"]}

                if new_attr["gene_id"] == new_attr["transcript_id"]:
                    new_attr["transcript_id"] = new_attr["transcript_id"] + "-T1"
                    new_attr["Parent"] = new_attr["Parent"] + "-T1"

                new_line = format_gff_line(fields, new_attr, feature)
                outfile.write(new_line)
            elif feature == "ncRNA":
                outfile.write(line)
            elif feature == "CDS":
                new_attr = {k: v for k, v in attr_dict.items() if k in ["ID", "Parent", "gene_id", "transcript_id", "exon_id"]}

                if new_attr["gene_id"] == new_attr["transcript_id"]:
                    new_attr["transcript_id"] = new_attr["transcript_id"] + "-T1"
                    new_attr["Parent"] = new_attr["Parent"] + "-T1"

                new_line = format_gff_line(fields, new_attr, feature)
                outfile.write(new_line)

            else:
                # Other features are ignored
                continue

def parse_attributes(attr_str):
    """
    Parse the attributes string into a dictionary.
    """
    attr_pairs = [attr.split('=') for attr in attr_str.split(';') if '=' in attr]
    return {k: v for k, v in attr_pairs}

def format_gff_line(fields, attributes, feature=None):
    """
    Format the fields and attributes into a GFF line.
    """
    attr_str = ';'.join(f'{k}={v}' for k, v in attributes.items())
    if feature:
        fields[2] = feature
    return '\t'.join(fields[:8] + [attr_str]) + '\n'

def modify_id(original_id, feature_type):
    """
    Modify the ID to a new format based on feature type.
    """
    # Example modification: append "_mod" and feature type to the original ID
    return f"{original_id}_mod_{feature_type}"

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python modify_gff.py <input_gff> <output_gff>")
        sys.exit(1)

    input_gff = sys.argv[1]
    output_gff = sys.argv[2]

    modify_gff(input_gff, output_gff)
