
import pandas as pd
import sys

def mod_ids(row, eq):
    attributes = row['attributes']
    if 'ID=' in attributes:
        id_b = attributes.split('ID=')[1].split(';')[0]
        if id_b in eq:
            id_a = eq[id_b]
            attributes = attributes.replace(f'ID={id_b}', f'ID={id_a}')
    if 'Parent=' in attributes:
        parent_b = attributes.split('Parent=')[1].split(';')[0]
        if parent_b in eq:
            parent_a = eq[parent_b]
            attributes = attributes.replace(f'Parent={parent_b}', f'Parent={parent_a}')
    return attributes

def filter_rows(row, ids):
    id_found = any(f'ID={id_str}' in row['attributes'] for id_str in ids)
    parent_found = any(f'Parent={id_str}' in row['attributes'] for id_str in ids)
    return not (id_found or parent_found)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python adjust_gff.py original.gff funannotate.gff output.gff")
        sys.exit(1)
    
    # Parsear ambos archivos
    
    original = pd.read_csv(sys.argv[1], sep='\t', comment='#', header=None, names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    funannotate = pd.read_csv(sys.argv[2], sep='\t', comment='#', header=None, names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

    output_gff = sys.argv[3]
    
    genes_a = original[original['type'] == 'gene']
    genes_b = funannotate[funannotate['type'] == 'gene']

    # If there is overlapping of more than one gene remove it from funannotate annotation

    eqGenes = {}
    genes_b_to_discard = []
    for _, gene_b in genes_b.iterrows():
        subset_a = genes_a[(genes_a['seqid'] == gene_b['seqid']) &
                           (genes_a['strand'] == gene_b['strand']) &
                           (genes_a['start'] >= gene_b['start']) &
                           (genes_a['end'] <= gene_b['end'])]
        if not subset_a.empty:
            # Si hay más de un gen en A que cumple la condición, nos quedamos con el primero
            if len(subset_a) > 1:
                genes_b_to_discard.append(gene_b['attributes'].split('ID=')[1].split(';')[0])
            else:
                gene_a = subset_a.iloc[0]
                id_b = gene_b['attributes'].split('ID=')[1].split(';')[0]
                id_a = gene_a['attributes'].split('ID=')[1].split(';')[0]
                eqGenes[id_b] = id_a

    # Now add the transcripts
    transcripts_b = funannotate[funannotate['type'] == 'mRNA']

    for _, transcript_b in transcripts_b.iterrows():
        id_b = transcript_b['attributes'].split('ID=')[1].split(';')[0]
        parent_b = transcript_b['attributes'].split('Parent=')[1].split(';')[0]

        if parent_b in eqGenes:
            transcriptNumber = id_b.split("-T")[1]
            new_id_b = eqGenes[parent_b] + "-T" + transcriptNumber
            eqGenes[id_b] = new_id_b
        if parent_b in genes_b_to_discard:
            genes_b_to_discard.append(id_b)


    # Now add the rest of the features

    other_bs = funannotate[funannotate['type'].isin(['five_prime_UTR', 'exon', 'three_prime_UTR', 'CDS'])]

    for _, other_b in other_bs.iterrows():
        id_b = other_b['attributes'].split('ID=')[1].split(';')[0]
        parent_b = other_b['attributes'].split('Parent=')[1].split(';')[0]

        if parent_b in eqGenes and other_b['type'] == 'exon':
            assignedNumber = id_b.split(".exon")[1]
            new_id_b = eqGenes[parent_b] + ".exon" + assignedNumber
            eqGenes[id_b] = new_id_b
        if parent_b in eqGenes and other_b['type'] == 'CDS':
            assignedNumber = id_b.split(".cds")[1]
            new_id_b = eqGenes[parent_b] + ".cds" + assignedNumber
            eqGenes[id_b] = new_id_b
        if parent_b in eqGenes and other_b['type'] == 'five_prime_UTR':
            assignedNumber = id_b.split(".utr5p")[1]
            new_id_b = eqGenes[parent_b] + ".utr5p" + assignedNumber
            eqGenes[id_b] = new_id_b
        if parent_b in eqGenes and other_b['type'] == 'three_prime_UTR':
            assignedNumber = id_b.split(".utr3p")[1]
            new_id_b = eqGenes[parent_b] + ".utr3p" + assignedNumber
            eqGenes[id_b] = new_id_b

    funannotate['attributes'] = funannotate.apply(mod_ids, axis=1, eq=eqGenes)
    filteredFunannotate = funannotate[funannotate.apply(filter_rows, ids=genes_b_to_discard, axis=1)]
    
    ####
    genes_a_to_add = original[~original['attributes'].str.contains('|'.join(eqGenes.values()))]
    final_funannotate = pd.concat([filteredFunannotate, genes_a_to_add], ignore_index=True)
    ####

    final_funannotate.to_csv(output_gff, sep='\t', index=False, header=False)
