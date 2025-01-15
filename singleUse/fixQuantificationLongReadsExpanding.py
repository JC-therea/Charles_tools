
import sys

def parse_attributes(attributes):
    """Parse the attributes column of a GTF line into a dictionary."""
    attributes_dict = {}
    for attr in attributes.split(';'):
        if attr.strip():
            key, value = attr.strip().split(' ')
            attributes_dict[key] = value.strip('"')
    return attributes_dict

def format_attributes(attributes_dict):
    """Format the attributes dictionary back into a GTF attributes string."""
    return '; '.join(f'{key} "{value}"' for key, value in attributes_dict.items()) + ';'

def process_gtf(input_gtf, output_gtf):
    reference_gene = {}
    reference_gene_expanded = {}
    reference_transcript = {}

    # Do a first loop to get all starting and end positions of the genes
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            attributes_dict = parse_attributes(attributes)

            if feature == 'gene':
                # Save the reference transcript's start and end positions
                gene_id = attributes_dict['gene_id']
                reference_gene_expanded[gene_id] = (int(start), int(end))

    # Now read all the reference transcripts
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            attributes_dict = parse_attributes(attributes)

            if feature == 'transcript' and source == 'Coding':
                
                gene_id = attributes_dict['gene_id']
                transcript_id = attributes_dict['transcript_id']
                
                if gene_id in reference_gene_expanded:
                    
                    reference_transcript[transcript_id] = (int(start), int(end))
                    ref_start, ref_end = reference_gene_expanded[gene_id]
                    new_start = min(int(start), ref_start)
                    new_end = max(int(end), ref_end)
                    
                    reference_gene_expanded[gene_id] = (new_start, new_end)                        
    # Now write the output file
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            attributes_dict = parse_attributes(attributes)

            if feature == 'transcript' and source == 'IsoQuant':
                outfile.write(line)
            elif feature == 'transcript' and source == 'Coding':

                gene_id = attributes_dict['gene_id']

                if gene_id in reference_gene_expanded:
                    fields[3] = str(reference_gene_expanded[gene_id][0])
                    fields[4] = str(reference_gene_expanded[gene_id][1])
                    outfile.write('\t'.join(fields) + '\n')
            elif feature == 'exon' and source == 'Coding':
                gene_id = attributes_dict['gene_id']
                transcript_id = attributes_dict['transcript_id']

                # extend 5 or 3 ends if the gene is bigger and the reference was equal to the previous

                if int(start) == reference_transcript[transcript_id][0]:
                    if reference_gene_expanded[gene_id][0] < int(start):
                        ref_start = reference_gene_expanded[gene_id][0]
                        new_start = min(int(start), ref_start)
                        fields[3] = str(new_start)
                if int(end) == reference_transcript[transcript_id][1]:
                    if reference_gene_expanded[gene_id][1] > int(end):
                        ref_end = reference_gene_expanded[gene_id][1]
                        new_end = max(int(end), ref_end)
                        fields[4] = str(new_end)
                
                outfile.write('\t'.join(fields) + '\n')
            else:
                outfile.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python adjust_gtf.py input.gtf output.gtf")
        sys.exit(1)

    input_gtf = sys.argv[1]
    output_gtf = sys.argv[2]

    process_gtf(input_gtf, output_gtf)
