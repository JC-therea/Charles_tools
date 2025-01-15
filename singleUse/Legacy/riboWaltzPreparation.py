from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

def extend_region(start, end, seq_len, extension=50):
    """
    Extend a given region by a specified extension.
    """
    new_start = max(start - extension, 1)
    new_end = min(end + extension, seq_len)
    return new_start, new_end

def create_mRNA_exons(gene, seq_dict, extension=50):
    """
    Given a gene, create mRNA and exon features for the CDS with
    a specified extension.
    """
    seq_id = gene.seqid
    seq_len = len(seq_dict[seq_id])
    strand = gene.strand
    gene_start = gene.start
    gene_end = gene.end
    gene_id = gene.id
    gene_name = gene.attributes['Name'][0]
    gene_mRNA_id = "{}.mRNA".format(gene_id)
    gene_exon_id = "{}.exon".format(gene_id)

    # Extend the CDS region by the specified extension
    start, end = extend_region(gene_start, gene_end, seq_len, extension)

    # Get the CDS region
    cds_feature = None
    for feature in gene.sub_features:
        if feature.type == "CDS":
            cds_feature = feature
            break

    if cds_feature:
        cds_start = cds_feature.start
        cds_end = cds_feature.end
        cds_seq = str(cds_feature.extract(seq_dict[seq_id].seq))

        # Create the mRNA feature
        mRNA_feature = GFF.FeatureLocation(start, end, strand)
        mRNA_attributes = {"ID": [gene_mRNA_id], "Name": [gene_name]}
        mRNA_feature = GFF.feature_from_seqrecord(SeqRecord(Seq(""), id=seq_id), "mRNA", mRNA_feature, attributes=mRNA_attributes)

        # Create the exon feature
        exon_start = start if strand == 1 else end
        exon_end = cds_start if strand == 1 else cds_end
        exon_feature = GFF.FeatureLocation(exon_start, exon_end, strand)
        exon_attributes = {"ID": [gene_exon_id]}
        exon_feature = GFF.feature_from_seqrecord(SeqRecord(Seq(""), id=seq_id), "exon", exon_feature, attributes=exon_attributes)

        # Add the UTR regions
        if cds_start > start:
            utr_start = start if strand == 1 else cds_end + 1
            utr_end = cds_start - 1 if strand == 1 else end
            utr_feature = GFF.FeatureLocation(utr_start, utr_end, strand)
            utr_attributes = {"ID": ["{}-UTR5".format(gene_id)]}
            utr_feature = GFF.feature_from_seqrecord(SeqRecord(Seq(""), id=seq_id), "UTR5", utr_feature, attributes=utr_attributes)
            gene.sub_features.append(utr_feature)

        if cds_end < end:
            utr_start = cds_end + 1 if strand == 1 else start
            utr_end = end if strand == 1 else cds_start - 1
            utr_feature = sda
# Add UTR regions to mRNAs
for mrna in mRNAs:
    gene_id = mrna['Parent']
    gene = genes[gene_id]
    for exon in mrna['exons']:
        if exon['strand'] == '+':
            if exon['start'] > gene['start']:
                utr5_start = gene['start']
                utr5_end = exon['start'] - 1
                mrna['UTR5'].append({'start': utr5_start, 'end': utr5_end, 'strand': '+'})
        else:
            if exon['end'] < gene['end']:
                utr5_start = exon['end'] + 1
                utr5_end = gene['end']
                mrna['UTR5'].append({'start': utr5_start, 'end': utr5_end, 'strand': '-'})
    mrna['UTR5'] = sorted(mrna['UTR5'], key=lambda x: x['start'])
    mrna['UTR3'] = []
    for i in range(len(mrna['exons']) - 1):
        exon = mrna['exons'][i]
        next_exon = mrna['exons'][i+1]
        if exon['strand'] == '+':
            utr3_start = exon['end'] + 1
            utr3_end = next_exon['start'] - 1
            mrna['UTR3'].append({'start': utr3_start, 'end': utr3_end, 'strand': '+'})
        else:
            utr3_start = next_exon['end'] + 1
            utr3_end = exon['start'] - 1
            mrna['UTR3'].append({'start': utr3_start, 'end': utr3_end, 'strand': '-'})
    if mrna['strand'] == '+':
        last_exon = mrna['exons'][-1]
        if last_exon['end'] < gene['end']:
            utr3_start = last_exon['end'] + 1
            utr3_end = min(gene['end'], last_exon['end'] + 50)
            mrna['UTR3'].append({'start': utr3_start, 'end': utr3_end, 'strand': '+'})
    else:
        first_exon = mrna['exons'][0]
        if first_exon['start'] > gene['start']:
            utr3_start = max(gene['start'], first_exon['start'] - 50)
            utr3_end = first_exon['start'] - 1
            mrna['UTR3'].append({'start': utr3_start, 'end': utr3_end, 'strand': '-'})
    mrna['UTR3'] = sorted(mrna['UTR3'], key=lambda x: x['start'])

# Write output to file
with open(output_file, 'w') as f:
    for gene_id in genes:
        gene = genes[gene_id]
        print_gff(gene, 'gene', f)
        for mrna in gene['mRNAs']:
            print_gff(mrna, 'mRNA', f)
            for exon in mrna['exons']:
                print_gff(exon, 'exon', f)
            for utr in mrna['UTR5']:
                print_gff(utr, 'UTR5', f)
            for utr in mrna['UTR3']:
                print_gff


