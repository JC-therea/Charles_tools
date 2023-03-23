from Bio import SeqIO
import argparse


def get_args():
    parser = argparse.ArgumentParser(description='Extract features from a genome in FASTA format based on a GFF file')
    parser.add_argument('-gff', '--gff_file', required=True, help='Path to the GFF file')
    parser.add_argument('-fasta', '--fasta_file', required=True, help='Path to the genome FASTA file')
    parser.add_argument('-f', '--features', nargs='+', required=True, help='Features to extract')
    parser.add_argument('-o', '--output_file', required=True, help='Output file in FASTA format')
    return parser.parse_args()


def get_feature_seqs(gff_file, fasta_file, features, output_file):
    # parse the GFF file
    feature_seqs = {}
    with open(gff_file, 'r') as gff:
        for line in gff:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if fields[2] in features:
                    # If its a UTR I want to get not the ID but the parent name
                    if "UTR" in fields[2]:
                        feature_id = fields[8].split(';')[1].split('=')[1]
                    else:
                        feature_id = fields[8].split(';')[0].split('=')[1]
                    if feature_id not in feature_seqs:
                        feature_seqs[feature_id] = ''
                    seq_start = int(fields[3]) - 1  # GFF is 1-indexed, but Python is 0-indexed
                    seq_end = int(fields[4])
                    seq_strand = fields[6]
                    with open(fasta_file, 'r') as genome:
                        for record in SeqIO.parse(genome, 'fasta'):
                            if record.id == fields[0]:
                                if seq_strand == '-':
                                    feature_seq = str(record.seq[seq_start:seq_end].reverse_complement())
                                else:
                                    feature_seq = str(record.seq[seq_start:seq_end])
                                feature_seqs[feature_id] += feature_seq

    # write the output file
    with open(output_file, 'w') as out:
        for feature_id, seq in feature_seqs.items():
            out.write(f'>{feature_id}\n{seq}\n')


if __name__ == '__main__':
    args = get_args()
    get_feature_seqs(args.gff_file, args.fasta_file, args.features, args.output_file)
