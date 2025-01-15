import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Report the number of times each singlet appears in a given multi fasta file')
    parser.add_argument('-fasta', '--fasta_file', required=True, help='Path to the multi fasta file')
    parser.add_argument('-o', '--output_file', required=True, help='Output file in FASTA format')
    return parser.parse_args()

def count_singlets(sequence):
    singlets = ['A','C','T','G']
    singlet_counts = {singlet: 0 for singlet in singlets}

    for i in range(len(sequence)):
        singlet = sequence[i]
        if singlet in singlet_counts:
            singlet_counts[singlet] += 1

    return singlet_counts


def analyze_singlets_in_multifasta(multifasta_file, output_file):
    with open(output_file, 'w') as output:
        output.write("Gene\tNucleotide\tCount\tSequence_Length\n")
        with open(multifasta_file, 'r') as fasta_file:
            gene = ""
            sequence = ""
            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    if gene and sequence:
                        singlet_counts = count_singlets(sequence)
                        seq_length = len(sequence)
                        for singlet, count in singlet_counts.items():
                            output.write(f"{gene}\t{singlet}\t{count}\t{seq_length}\n")
                    gene = line[1:]
                    sequence = ""
                else:
                    sequence += line.upper()
            if gene and sequence:
                singlet_counts = count_singlets(sequence)
                seq_length = len(sequence)
                for singlet, count in singlet_counts.items():
                    output.write(f"{gene}\t{singlet}\t{count}\t{seq_length}\n")


args = get_args()
analyze_singlets_in_multifasta(args.fasta_file, args.output_file)
