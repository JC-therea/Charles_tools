

# Función para leer el bedgraph en un DataFrame
def read_coverage(file):
    return pd.read_csv(file, sep='\t', header=None, names=['chrom', 'pos', 'coverage'])

# Función para leer el GFF y extraer las coordenadas de las CDS
def gff2pd(file):
    regions = []
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                chrom = parts[0]
                feature = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attr = parts[8]
                if feature != 'CDS':
                    gene_id = [field.split("=")[1] for field in attr.split(";") if field.startswith("ID")][0]
                if feature == 'CDS':
                    gene_id = attr.split("Parent=")[1].split(";")[0]
                regions.append({"chrom": chrom, "feature": feature, "start": start, "end": end, "strand": strand, "gene_id": gene_id})
    return pd.DataFrame(regions)

def geneId2transcriptId(file):
    geneId2transcriptId = {}
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                chrom = parts[0]
                feature = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attr = parts[8]
                id = [field.split("=")[1] for field in attr.split(";") if field.startswith("ID")][0]
                if feature == "mRNA" or feature == "transcript":
                    gene_id = attr.split("Parent=")[1].split(";")[0]
                    geneId2transcriptId[gene_id] = id
    return geneId2transcriptId


# Función para calcular el pico máximo en cada CDS
def calculate_max_coverage(gff_df, coverage_df):
    if gff_df['feature'].str.contains('transcript', regex=True).any():
        transcriptName = "transcript"
    if gff_df['feature'].str.contains('mRNA', regex=True).any():
        transcriptName = "mRNA"

    transcripts_df = gff_df[gff_df['feature'] == transcriptName]
    coverage_per_transcript = []

    start_df = transcripts_df[['chrom', 'feature', 'start', 'gene_id', 'strand']].merge(
        coverage_df, left_on=['chrom', 'start'], right_on=['chrom', 'pos'], how='left'
    ).rename(columns={'coverage': 'start_coverage'})

    end_df = transcripts_df[['chrom', 'feature', 'end', 'gene_id', 'strand']].merge(
        coverage_df, left_on=['chrom', 'end'], right_on=['chrom', 'pos'], how='left'
    ).rename(columns={'coverage': 'end_coverage'})

    merged_df = pd.merge(start_df[['gene_id', 'feature', 'chrom', 'start', 'strand', 'start_coverage']],
                         end_df[['gene_id', 'end', 'end_coverage']], on='gene_id', how='left')
    return merged_df

def extend_transcripts(gff_df_coverage, coverage_df, threshold=0.1, min_coverage=5):
    extended_transcripts = []

    # Agrupar el archivo de cobertura por cromosoma para reducir el filtrado en cada iteración
    grouped_coverage = coverage_df.groupby('chrom')

    for index, row in gff_df_coverage.iterrows():
        chrom, start, start_coverage, end, end_coverage = row['chrom'], row['start'], row['start_coverage'], row['end'], row['end_coverage']

        # Determinar umbral de cobertura
        # New part here!
        start_threshold = max(start_coverage * threshold, min_coverage)
        start_upper_treshold = start_coverage * 2
        end_threshold = max(end_coverage * threshold, min_coverage)
        end_upper_treshold = end_coverage * 2

        # Filtrar datos de cobertura una sola vez por cada cromosoma
        if chrom in grouped_coverage.groups:
            chrom_coverage_df = grouped_coverage.get_group(chrom)

            # Obtener posiciones de inicio y fin extendidas basadas en el umbral
            start_filtered = chrom_coverage_df.loc[((chrom_coverage_df['coverage'] < start_threshold) | (chrom_coverage_df['coverage'] > start_upper_treshold )) & (chrom_coverage_df['pos'] < start), 'pos']
            end_filtered = chrom_coverage_df.loc[((chrom_coverage_df['coverage'] < end_threshold) | (chrom_coverage_df['coverage'] > end_upper_treshold )) & (chrom_coverage_df['pos'] > end), 'pos']

            extended_start = start_filtered.max() if not start_filtered.empty else start
            extended_end = end_filtered.min() if not end_filtered.empty else end

            if extended_start == start - 1:
                extended_start = start
            if extended_end == end + 1:
                extended_end = end

            extended_transcripts.append({
                "gene_id": row['gene_id'],
                "feature": row['feature'],
                "chrom": chrom,
                "start": start,
                "end": end,
                "extended_start": extended_start,
                "extended_end": extended_end,
                "strand": row['strand'],
                "start_coverage": start_coverage,
                "end_coverage": end_coverage
            })

    return pd.DataFrame(extended_transcripts)

# Función para acortar genes, mRNA y exones en caso de solapamiento con CDS en misma hebra

def shorten_overlapping_features(gffExpanded,gff_df):
    # Filtrar solo las filas relevantes de CDS y otros features
    cds_df = gff_df[gff_df['feature'] == 'CDS'][['chrom', 'start', 'end', 'strand', 'gene_id']]
    grouped_cds = cds_df.groupby('chrom')
    #cds_df.to_csv("cds.tsv", sep="\t")

    shortenGff = []

    for index, row in gffExpanded.iterrows():
        chrom, start, extended_start, end, extended_end, strand, gene_id = row['chrom'], row['start'], row['extended_start'], row['end'], row['extended_end'], row['strand'], row['gene_id']
        updated_start = 0
        updated_end = 0

        if chrom in grouped_cds.groups:
            chrom_cds_df = grouped_cds.get_group(chrom)

            # Filtrar CDS que están antes del `extended_start` y después del `extended_end` en la misma hebra
            mask_5end = ((chrom_cds_df['end'] > extended_start) &
                (chrom_cds_df['end'] < start) &
                (chrom_cds_df['gene_id'] != gene_id) &
                (chrom_cds_df['strand'] == strand)
            )
            mask_3end = ((chrom_cds_df['start'] < extended_end) &
                (chrom_cds_df['start'] > end) &
                (chrom_cds_df['gene_id'] != gene_id) &
                (chrom_cds_df['strand'] == strand)
            )
            filtered_5end = chrom_cds_df['end'][mask_5end].max()
            filtered_3end = chrom_cds_df['start'][mask_3end].min()

            # Obtener el valor máximo de inicio y el mínimo de fin si hay resultados válidos
            if not pd.isna(filtered_5end) and extended_start != start:
                updated_start = filtered_5end + 1
            if not pd.isna(filtered_3end) and extended_end != end:
                updated_end = filtered_3end - 1

        final_start = extended_start
        final_end = extended_end
        if updated_start != 0:
            final_start = updated_start
        # Actualizar el DataFrame expandido con los nuevos valores de inicio y fin
        if updated_end != 0:
            final_end = updated_end

        shortenGff.append({
            "gene_id": row['gene_id'],
            "feature": row['feature'],
            "chrom": chrom,
            "start": start,
            "end": end,
            "extended_start": final_start,
            "extended_end": final_end,
            "strand": row['strand'],
        })

    return pd.DataFrame(shortenGff)

def remove_completly_overlapping_extensions(gffExpanded_defined,gff_df):
    # Filtrar solo las filas relevantes de CDS y otros features
    cds_df = gff_df[gff_df['feature'] == 'CDS'][['chrom', 'start', 'end', 'strand', 'gene_id']]
    grouped_cds = cds_df.groupby('chrom')
    #cds_df.to_csv("cds.tsv", sep="\t")

    finalGff = []

    for index, row in gffExpanded_defined.iterrows():
        chrom, start, extended_start, end, extended_end, strand, gene_id = row['chrom'], row['start'], row['extended_start'], row['end'], row['extended_end'], row['strand'], row['gene_id']


        if chrom in grouped_cds.groups:
            chrom_cds_df = grouped_cds.get_group(chrom)

            # If covers the start or end point of the CDS

            mask_5end = ((extended_start <= chrom_cds_df['start']) &
                (start > chrom_cds_df['start']) &
                (chrom_cds_df['gene_id'] != gene_id) &
                (chrom_cds_df['strand'] == strand)
            )

            mask_3end = ((extended_end >= chrom_cds_df['end']) &
                (end < chrom_cds_df['end']) &
                (chrom_cds_df['gene_id'] != gene_id) &
                (chrom_cds_df['strand'] == strand)
            )

            # If it's fully covered by another CDS

            mask_5end_fully = ((extended_start > chrom_cds_df['start']) &
                (extended_start < chrom_cds_df['end']) &
                (chrom_cds_df['gene_id'] != gene_id) &
                (chrom_cds_df['strand'] == strand)
            )

            mask_3end_fully = ((extended_end > chrom_cds_df['start']) &
                (extended_end < chrom_cds_df['end']) &
                (chrom_cds_df['gene_id'] != gene_id) &
                (chrom_cds_df['strand'] == strand)
            )

            if mask_5end.any() or mask_5end_fully.any():
                finalExtendedStart = start
            else:
                finalExtendedStart = extended_start

            if mask_3end.any() or mask_3end_fully.any():
                finalExtendedEnd = end
            else:
                finalExtendedEnd = extended_end

        finalGff.append({
            "gene_id": row['gene_id'],
            "feature": row['feature'],
            "chrom": chrom,
            "start": start,
            "end": end,
            "extended_start": finalExtendedStart,
            "extended_end": finalExtendedEnd,
            "strand": row['strand'],
        })

    return pd.DataFrame(finalGff)


def update_gff_from_extended(gff_file, extended_df,gene2Trans, output_file):
    if extended_df['feature'].str.contains('transcript', regex=True).any():
        transcriptName = "transcript"
    if extended_df['feature'].str.contains('mRNA', regex=True).any():
        transcriptName = "mRNA"
    with open(gff_file, 'r') as infile, open(output_file, 'w+') as outfile:
        for line in infile:

            # Escribir comentarios y encabezado sin modificar
            if line.startswith("#"):
                outfile.write(line)
                continue

            parts = line.strip().split("\t")
            chrom, feature, start, end, strand, attr = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6], parts[8]
            attributes = dict(field.split("=") for field in attr.split(";"))

            # Identificar el ID adecuado para buscar en `extended_df`
            if feature == "gene":
                if attributes["ID"] not in gene2Trans.keys():
                    continue
                gene_id = gene2Trans[attributes["ID"]]  # ID sin ".mrna"
            elif feature == transcriptName:
                gene_id = attributes["ID"]  # ID con ".mrna"
            elif feature == "exon":
                gene_id = attributes["Parent"]  # transcript_id para exones
            else:
                # Si no es un feature relevante, escribir la línea original
                outfile.write(line)
                continue

            # Buscar las posiciones extendidas en `extended_df`
            extended_row = extended_df[(extended_df["gene_id"] == gene_id) & (extended_df["chrom"] == chrom)]
            if not extended_row.empty:
                if extended_row["start"].values[0] == start:
                    new_start = extended_row["extended_start"].values[0]
                else:
                    new_start = start
                if extended_row["end"].values[0] == end:
                    new_end = extended_row["extended_end"].values[0]
                else:
                    new_end = end
                parts[3], parts[4] = str(new_start), str(new_end)

            # Escribir la línea actualizada en el archivo de salida
            outfile.write("\t".join(parts) + "\n")

if __name__ == '__main__':

    from datetime import datetime
    import pandas as pd
    import sys

    try:
        functionName, covPlusPath, covMinusPath, gff, threshold, minReads, outPath  = sys.argv
        threshold = float(threshold)
        minReads = int(minReads)
    except:
        print("error: utr.bam.extender.py covPlusPath, covMinusPath, gff, threshold of coverage (in %), minReads, outPath")
        sys.exit()

    # Read files
    print("Reading positive bedgraph...")
    start_time = datetime.now()
    positive_bedgraph_df = read_coverage(covPlusPath)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    start_time = datetime.now()
    print("Reading negative bedgraph...")
    negative_bedgraph_df = read_coverage(covMinusPath)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    start_time = datetime.now()
    print("Reading gff file...")
    gff_df = gff2pd(gff)
    g2t = geneId2transcriptId(gff)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

    # Estimating peaks at each side

    start_time = datetime.now()
    print("Estimating peaks at 5' and 3' per transcript + strand...")
    positive_max_coverage_df = calculate_max_coverage(gff_df[gff_df['strand'] == '+'], positive_bedgraph_df)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    start_time = datetime.now()
    print("Estimating peaks at 5' and 3' per transcript - strand...")
    negative_max_coverage_df = calculate_max_coverage(gff_df[gff_df['strand'] == '-'], negative_bedgraph_df)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

    # Extending transcripts until first position under the threshold

    start_time = datetime.now()
    print("Extending transcripts according to coverage + strand...")
    positive_extended_df = extend_transcripts(positive_max_coverage_df, positive_bedgraph_df, threshold = threshold, min_coverage = minReads)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    start_time = datetime.now()
    print("Extending transcripts according to coverage - strand...")
    negative_extended_df = extend_transcripts(negative_max_coverage_df, negative_bedgraph_df, threshold = threshold, min_coverage = minReads)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    # Concatenate both strands
    extended_df = pd.concat([positive_extended_df, negative_extended_df])

    # Limit the extension to other CDSs

    start_time = datetime.now()
    print("Adjusting UTRs for not overlapping CDSs...")
    adjusted_df = shorten_overlapping_features(extended_df, gff_df)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

    # Remove UTRs that still overlapping CDSs

    start_time = datetime.now()
    print("Removing remaining UTRs overlapping CDSs...")
    final_df = remove_completly_overlapping_extensions(adjusted_df, gff_df)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))


    # Translate from pandas dataframe to gff

    start_time = datetime.now()
    print("Updating gff...")
    update_gff_from_extended(gff, final_df, g2t, outPath)
    end_time = datetime.now()
    print('Done in: {}'.format(end_time - start_time))
