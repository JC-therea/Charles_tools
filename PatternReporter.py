import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(description='Report the number of times a pattern appears in a given multi fasta file')
    parser.add_argument('-p', '--pattern', required=True, help='Pattern to search in muti fasta file')
    parser.add_argument('-fasta', '--fasta_file', required=True, help='Path to the multi fasta file')
    parser.add_argument('-o', '--output_file', required=True, help='Output file in FASTA format')
    return parser.parse_args()


def search_pattern(multifasta_file, output_file, pattern):
    sequences = {}  # Dictionary to store sequence names and sequences

    # Read multifasta file and store sequences
    with open(multifasta_file, 'r') as fasta:
        current_seq = ''
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                current_seq = line[1:]
                sequences[current_seq] = ''
            else:
                sequences[current_seq] += line.upper()

    pattern_count = {}  # Dictionary to store pattern counts for each sequence

    # Initialize pattern count for each sequence
    for seq_name in sequences:
        pattern_count[seq_name] = 0

    # Search for pattern in each sequence
    for seq_name, sequence in sequences.items():
        # Convert 'N' to a regex pattern that matches any letter (A/C/T/G)
        seq_pattern = pattern.replace('N', '[ACGT]')

        # Count the number of non-overlapping pattern occurrences in the sequence
        count = len(re.findall(f'(?=({seq_pattern}))', sequence))
        pattern_count[seq_name] = count

    # Write pattern counts to the output file
    with open(output_file, 'w') as output:
        output.write("Gene\tPattern count\n")
        for seq_name, count in pattern_count.items():
            output.write(f"{seq_name}\t{count}\n")

    print("Pattern search completed. Results written to", output_file)


args = get_args()
search_pattern(args.fasta_file, args.output_file, args.pattern)
