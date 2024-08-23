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

def process_gtf(input_gtf, reference_gtf, output_gtf):
    reference_transcripts = {}
    isoform_transcripts = {}
    equivalence_ref = {}

    # Do a first loop to get all starting and end positions
    with open(reference_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            attributes_dict = parse_attributes(attributes)

            if feature == 'transcript':
                # Save the reference transcript's start and end positions
                transcript_id = attributes_dict['transcript_id']
                reference_transcripts[transcript_id] = (int(start), int(end))

    # Now read all the isoQuant transcripts
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            attributes_dict = parse_attributes(attributes)

            if feature == 'transcript' and source == 'IsoQuant':
                if 'similar_reference_id' in attributes_dict:

                    transcript_id = attributes_dict['transcript_id']
                    reference_id = attributes_dict['similar_reference_id']

                    if reference_id in reference_transcripts:
                        equivalence_ref[transcript_id] = reference_id

                        ref_start, ref_end = reference_transcripts[reference_id]

                        # Adjust the start and end positions if necessary
                        new_start = max(int(start), ref_start)
                        new_end = min(int(end), ref_end)

                        isoform_transcripts[transcript_id] = (new_start, new_end)

    # Now write the output file
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            attributes_dict = parse_attributes(attributes)

            if feature == 'transcript' and source == 'Coding':
                outfile.write(line)
            elif feature == 'transcript' and source == 'IsoQuant':
                if 'similar_reference_id' in attributes_dict.keys():

                    transcript_id = attributes_dict['transcript_id']
                    reference_id = attributes_dict['similar_reference_id']

                    if reference_id in reference_transcripts:

                        fields[3] = str(isoform_transcripts[transcript_id][0])
                        fields[4] = str(isoform_transcripts[transcript_id][1])
                        outfile.write('\t'.join(fields) + '\n')
            elif feature == 'gene':

                reference_id = attributes_dict['gene_id'].replace("-G1", "")

                if reference_id in reference_transcripts:

                    fields[3] = str(reference_transcripts[reference_id][0])
                    fields[4] = str(reference_transcripts[reference_id][1])
                    outfile.write('\t'.join(fields) + '\n')
            elif feature == 'exon' and source == 'IsoQuant':
                transcript_id = attributes_dict['transcript_id']

                if transcript_id in isoform_transcripts:
                    ref_start = isoform_transcripts[transcript_id][0]
                    ref_end = isoform_transcripts[transcript_id][1]
                    if int(end) < ref_start:
                        continue
                    elif int(start) > ref_end:
                        continue
                    else:
                        # Adjust the start and end positions if necessary
                        new_start = max(int(start), ref_start)
                        new_end = min(int(end), ref_end)

                        fields[3] = str(new_start)
                        fields[4] = str(new_end)
                        outfile.write('\t'.join(fields) + '\n')
            else:
                outfile.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python adjust_gtf.py input.gtf reference.gtf output.gtf")
        sys.exit(1)

    input_gtf = sys.argv[1]
    reference_gtf = sys.argv[2]
    output_gtf = sys.argv[3]

    process_gtf(input_gtf, reference_gtf,  output_gtf)
