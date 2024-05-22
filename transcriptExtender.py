import pandas as pd
import argparse

def read_gtf(gtf_file):
    """Reads a GTF file into a pandas DataFrame."""
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', names=columns)
    return gtf

def extend_features(gtf, extend_5utr, extend_3utr):
    """Extends 5' and 3' UTRs by specified sizes if they do not exist, and adjusts exons, genes, and transcripts."""
    features = ['gene', 'transcript', 'exon']
    transcriptOriCoords = {}
    for feature in features:
        for idx, row in gtf[gtf['feature'] == feature].iterrows():
            if feature == 'transcript':
                transcript = row['attribute'].split('transcript_id "')[1].split('";')[0]
                transcriptOriCoords[transcript] = [row['start'], row['end']]
            if feature != 'exon':
                if row['strand'] == '+':
                    gtf.at[idx, 'start'] = max(1, row['start'] - extend_5utr)
                    gtf.at[idx, 'end'] = row['end'] + extend_3utr
                else:
                    gtf.at[idx, 'end'] = row['end'] + extend_5utr
                    gtf.at[idx, 'start'] = max(1, row['start'] - extend_3utr)
            else:
                transcript = row['attribute'].split('transcript_id "')[1].split('";')[0]
                if row['strand'] == '+':
                    if row['start'] == transcriptOriCoords[transcript][0]:
                        gtf.at[idx, 'start'] = max(1, row['start'] - extend_5utr)
                    if row['end'] == transcriptOriCoords[transcript][1]:
                        gtf.at[idx, 'end'] = row['end'] + extend_3utr
                else:
                    if row['start'] == transcriptOriCoords[transcript][0]:
                        gtf.at[idx, 'start'] = max(1, row['start'] - extend_3utr)
                    if row['end'] == transcriptOriCoords[transcript][1]:
                        gtf.at[idx, 'end'] = row['end'] + extend_5utr

    return gtf

def extend_utr(gtf, extend_5utr, extend_3utr):
    """Extends 5' and 3' UTRs by specified sizes if they do not exist."""
    cds = gtf[gtf['feature'] == 'CDS']
    
    for i, row in cds.iterrows():
        if row['strand'] == '+':
            if gtf[(gtf['feature'] == 'five_prime_utr') & (gtf['attribute'] == row['attribute'])].empty:
                new_start = max(1, row['start'] - extend_5utr)
                new_row = row.copy()
                new_row['feature'] = 'five_prime_utr'
                new_row['start'] = new_start
                new_row['end'] = row['start'] - 1
                gtf = pd.concat([gtf,new_row], ignore_index=True).copy()
            if gtf[(gtf['feature'] == 'three_prime_utr') & (gtf['attribute'] == row['attribute'])].empty:
                new_end = row['end'] + extend_3utr
                new_row = row.copy()
                new_row['feature'] = 'three_prime_utr'
                new_row['start'] = row['end'] + 1
                new_row['end'] = new_end
                gtf = pd.concat([gtf,new_row], ignore_index=True).copy()
        else:
            if gtf[(gtf['feature'] == 'five_prime_utr') & (gtf['attribute'] == row['attribute'])].empty:
                new_end = row['end'] + extend_5utr
                new_row = row.copy()
                new_row['feature'] = 'five_prime_utr'
                new_row['start'] = row['end'] + 1
                new_row['end'] = new_end
                gtf = pd.concat([gtf,new_row], ignore_index=True).copy()
            if gtf[(gtf['feature'] == 'three_prime_utr') & (gtf['attribute'] == row['attribute'])].empty:
                new_start = max(1, row['start'] - extend_3utr)
                new_row = row.copy()
                new_row['feature'] = 'three_prime_utr'
                new_row['start'] = new_start
                new_row['end'] = row['start'] - 1
                gtf = pd.concat([gtf,new_row], ignore_index=True).copy()
    
    gtf = gtf.sort_values(by=['seqname', 'start', 'end']).reset_index(drop=True)
    return gtf

def write_gtf(gtf, output_file):
    """Writes the GTF DataFrame to a file."""
    gtf.to_csv(output_file, sep='\t', index=False, header=False, quoting=3)

def main():
    parser = argparse.ArgumentParser(description="Extend 5' and 3' UTRs in a GTF file.")
    parser.add_argument("gtf_file", help="Input GTF file")
    parser.add_argument("--extend_5utr", type=int, default=50, help="Size to extend 5' UTR (default: 50)")
    parser.add_argument("--extend_3utr", type=int, default=50, help="Size to extend 3' UTR (default: 50)")
    parser.add_argument("--output_file", default="extended_utrs.gtf", help="Output GTF file (default: extended_utrs.gtf)")
    
    args = parser.parse_args()
    
    gtf = read_gtf(args.gtf_file)
    gtf = extend_features(gtf, args.extend_5utr, args.extend_3utr)
    #extended_gtf = extend_utr(gtf, args.extend_5utr, args.extend_3utr)
    write_gtf(gtf, args.output_file)
    print(f"Extended UTRs written to {args.output_file}")

if __name__ == "__main__":
    main()

