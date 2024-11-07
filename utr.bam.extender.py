

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
                gene_id = [field.split("=")[1] for field in attr.split(";") if field.startswith("ID")][0]
                regions.append({"chrom": chrom, "feature": feature, "start": start, "end": end, "strand": strand, "gene_id": gene_id})
    return pd.DataFrame(regions)

# Función para calcular el pico máximo en cada CDS
def calculate_max_coverage(gff_df, coverage_df):
    coverage_per_transcript = []
    transcripts_df = gff_df[gff_df['feature'] == 'transcript']

    start_df = transcripts_df[['chrom', 'start', 'gene_id', 'strand']].merge(
        coverage_df, left_on=['chrom', 'start'], right_on=['chrom', 'pos'], how='left'
    ).rename(columns={'coverage': 'start_coverage'})

    end_df = transcripts_df[['chrom', 'end', 'gene_id', 'strand']].merge(
        coverage_df, left_on=['chrom', 'end'], right_on=['chrom', 'pos'], how='left'
    ).rename(columns={'coverage': 'end_coverage'})

    merged_df = pd.merge(start_df[['gene_id', 'chrom', 'start', 'strand', 'start_coverage']],
                         end_df[['gene_id', 'end', 'end_coverage']], on='gene_id', how='left')
    return merged_df

def extend_transcripts(gff_df_coverage, coverage_df, threshold=0.1, min_coverage=5):
    extended_transcripts = []

    # Agrupar el archivo de cobertura por cromosoma para reducir el filtrado en cada iteración
    grouped_coverage = coverage_df.groupby('chrom')

    for index, row in gff_df_coverage.iterrows():
        chrom, start, start_coverage, end, end_coverage = row['chrom'], row['start'], row['start_coverage'], row['end'], row['end_coverage']

        # Determinar umbral de cobertura
        start_threshold = max(start_coverage * threshold, min_coverage)
        end_threshold = max(end_coverage * threshold, min_coverage)

        # Filtrar datos de cobertura una sola vez por cada cromosoma
        if chrom in grouped_coverage.groups:
            chrom_coverage_df = grouped_coverage.get_group(chrom)

            # Obtener posiciones de inicio y fin extendidas basadas en el umbral
            start_filtered = chrom_coverage_df.loc[(chrom_coverage_df['coverage'] < start_threshold) & (chrom_coverage_df['pos'] < start), 'pos']
            end_filtered = chrom_coverage_df.loc[(chrom_coverage_df['coverage'] < end_threshold) & (chrom_coverage_df['pos'] > end), 'pos']

            extended_start = start_filtered.max() if not start_filtered.empty else start
            extended_end = end_filtered.min() if not end_filtered.empty else end

            if extended_start == start - 1:
                extended_start = start
            if extended_end == end + 1:
                extended_end = end

            extended_transcripts.append({
                "gene_id": row['gene_id'],
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

def update_gff_from_extended(gff_file, extended_df, output_file):
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
                gene_id = attributes["ID"] + ".mrna"  # ID sin ".mrna"
            elif feature == "transcript":
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
        print("error: utr.bam.extender.py covPlusPath, covMinusPath, gff, outPath")
        sys.exit()

    # Leer archivos
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
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    # Calcular los picos para cada CDS en función de la cadena
    start_time = datetime.now()
    print("Estimating maximum coverage per transcript + strand...")
    positive_max_coverage_df = calculate_max_coverage(gff_df[gff_df['strand'] == '+'], positive_bedgraph_df)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    start_time = datetime.now()
    print("Estimating maximum coverage per transcript - strand...")
    negative_max_coverage_df = calculate_max_coverage(gff_df[gff_df['strand'] == '-'], negative_bedgraph_df)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

    # Combinar los resultados

    start_time = datetime.now()
    print("Extending CDSs according to coverage + strand...")
    positive_extended_df = extend_transcripts(positive_max_coverage_df, positive_bedgraph_df, threshold = threshold, min_coverage = minReads)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    start_time = datetime.now()
    print("Extending CDSs according to coverage - strand...")
    negative_extended_df = extend_transcripts(negative_max_coverage_df, negative_bedgraph_df, threshold = threshold, min_coverage = minReads)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

    extended_df = pd.concat([positive_extended_df, negative_extended_df])

    start_time = datetime.now()
    print("Updating gff...")
    update_gff_from_extended(gff, extended_df, outPath)
    end_time = datetime.now()
    print('Done in: {}'.format(end_time - start_time))

    #extended_df.to_csv(outPath, sep = "\t", header=False, index=False)
    #print(max_coverage_df)
