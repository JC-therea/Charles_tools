import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Report the number of times each triplet appears in a given multi fasta file')
    parser.add_argument('-fasta', '--fasta_file', required=True, help='Path to the multi fasta file')
    parser.add_argument('-o', '--output_file', required=True, help='Output file in FASTA format')
    return parser.parse_args()


from itertools import product


def count_triplets(sequence):
    triplets = [''.join(triplet) for triplet in product('ACGT', repeat=3)]
    triplet_counts = {triplet: 0 for triplet in triplets}

    for i in range(len(sequence) - 2):
        triplet = sequence[i:i + 3]
        if triplet in triplet_counts:
            triplet_counts[triplet] += 1

    return triplet_counts


def analyze_triplets_in_multifasta(multifasta_file, output_file):
    with open(output_file, 'w') as output:
        output.write("Gene\tTriplet\tCount\tSequence_Length\n")
        with open(multifasta_file, 'r') as fasta_file:
            gene = ""
            sequence = ""
            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    if gene and sequence:
                        triplet_counts = count_triplets(sequence)
                        seq_length = len(sequence)
                        for triplet, count in triplet_counts.items():
                            output.write(f"{gene}\t{triplet}\t{count}\t{seq_length}\n")
                    gene = line[1:]
                    sequence = ""
                else:
                    sequence += line.upper()
            if gene and sequence:
                triplet_counts = count_triplets(sequence)
                seq_length = len(sequence)
                for triplet, count in triplet_counts.items():
                    output.write(f"{gene}\t{triplet}\t{count}\t{seq_length}\n")


args = get_args()
analyze_triplets_in_multifasta(args.fasta_file, args.output_file)