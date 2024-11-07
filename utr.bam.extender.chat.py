

# Función para leer el bedgraph en un DataFrame
def read_bedgraph(file):
    return pd.read_csv(file, sep='\t', header=None, names=['chrom', 'start', 'end', 'coverage'])

# Función para leer el GFF y extraer las coordenadas de las CDS
def read_gff(file):
    cds_regions = []
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                if parts[2] == "CDS":
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    attr = parts[8]
                    gene_id = [field.split("=")[1] for field in attr.split(";") if field.startswith("ID")][0]
                    cds_regions.append({"chrom": chrom, "start": start, "end": end, "strand": strand, "gene_id": gene_id})
    return pd.DataFrame(cds_regions)

# Función para calcular el pico máximo en cada CDS
def calculate_max_coverage(cds_df, bedgraph_df):
    max_coverage_per_cds = []
    for index, cds in cds_df.iterrows():
        # Filtrar por las posiciones de CDS en el bedgraph
        overlap = bedgraph_df[(bedgraph_df['chrom'] == cds['chrom']) &
                              (bedgraph_df['start'] < cds['end']) &
                              (bedgraph_df['end'] > cds['start'])]
        overlap = overlap.copy()
        # Ajustar los valores para que estén dentro de los límites del CDS
        overlap['overlap_start'] = overlap[['start', 'end']].max(axis=1)
        overlap['overlap_end'] = overlap[['start', 'end']].min(axis=1)
        overlap['overlap_start'] = overlap['overlap_start'].clip(lower=cds['start'])
        overlap['overlap_end'] = overlap['overlap_end'].clip(upper=cds['end'])

        # Obtener el valor máximo de cobertura
        max_coverage = overlap['coverage'].max()
        max_coverage_per_cds.append({
            "gene_id": cds['gene_id'],
            "chrom": cds['chrom'],
            "start": cds['start'],
            "end": cds['end'],
            "strand": cds['strand'],
            "max_coverage": max_coverage
        })
    return pd.DataFrame(max_coverage_per_cds)

# Función para extender la región CDS según el umbral de cobertura
def extend_cds_regions(cds_df, bedgraph_df, threshold=0.1):
    extended_regions = []

    for index, row in cds_df.iterrows():
        max_coverage = row['max_coverage']
        threshold_coverage = max_coverage * threshold
        start, end = row['start'], row['end']

        # Extensión en la dirección 5' y 3' basándose en la hebra
        if row['strand'] == '+':
            # Extender en dirección 5' (inicio)
            for pos in range(start - 1, bedgraph_df['start'].min() - 1, -1):
                overlap = bedgraph_df[(bedgraph_df['chrom'] == row['chrom']) & (bedgraph_df['start'] <= pos) & (bedgraph_df['end'] > pos)]
                if overlap.empty or overlap.iloc[0]['coverage'] < threshold_coverage:
                    break
                start -= 1

            # Extender en dirección 3' (final)
            for pos in range(end + 1, bedgraph_df['end'].max() + 1):
                overlap = bedgraph_df[(bedgraph_df['chrom'] == row['chrom']) & (bedgraph_df['start'] < pos) & (bedgraph_df['end'] >= pos)]
                if overlap.empty or overlap.iloc[0]['coverage'] < threshold_coverage:
                    break
                end += 1

        else:  # Si la hebra es negativa
            # Extender en dirección 5' (final)
            for pos in range(end + 1, bedgraph_df['end'].max() + 1):
                overlap = bedgraph_df[(bedgraph_df['chrom'] == row['chrom']) & (bedgraph_df['start'] < pos) & (bedgraph_df['end'] >= pos)]
                if overlap.empty or overlap.iloc[0]['coverage'] < threshold_coverage:
                    break
                end += 1

            # Extender en dirección 3' (inicio)
            for pos in range(start - 1, bedgraph_df['start'].min() - 1, -1):
                overlap = bedgraph_df[(bedgraph_df['chrom'] == row['chrom']) & (bedgraph_df['start'] <= pos) & (bedgraph_df['end'] > pos)]
                if overlap.empty or overlap.iloc[0]['coverage'] < threshold_coverage:
                    break
                start -= 1

        # Añadir la región extendida al resultado
        extended_regions.append({
            "chrom": row['chrom'],
            "source": "extended",
            "feature": "CDS",
            "start": start,
            "end": end,
            "score": ".",
            "strand": row['strand'],
            "frame": ".",
            "attributes": row['attributes']
        })

    return pd.DataFrame(extended_regions)

if __name__ == '__main__':

    from datetime import datetime
    import pandas as pd
    import sys

    try:
        functionName, covPlusPath, covMinusPath, gff, outPath  = sys.argv
    except:
        print("error: utr.bam.extender.py covPlusPath, covMinusPath, gff, outPath")
        sys.exit()

    # Leer archivos
    print("Reading positive bedgraph...")
    start_time = datetime.now()
    positive_bedgraph_df = read_bedgraph(covPlusPath)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    start_time = datetime.now()
    print("Reading negative bedgraph...")
    negative_bedgraph_df = read_bedgraph(covMinusPath)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    start_time = datetime.now()
    print("Reading gff file...")
    cds_df = read_gff(gff)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    # Calcular los picos para cada CDS en función de la cadena
    start_time = datetime.now()
    print("Estimating maximum coverage per transcript + strand...")
    positive_max_coverage_df = calculate_max_coverage(cds_df[cds_df['strand'] == '+'], positive_bedgraph_df)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    start_time = datetime.now()
    print("Estimating maximum coverage per transcript - strand...")
    negative_max_coverage_df = calculate_max_coverage(cds_df[cds_df['strand'] == '-'], negative_bedgraph_df)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

    # Combinar los resultados
    max_coverage_df = pd.concat([positive_max_coverage_df, negative_max_coverage_df])
    start_time = datetime.now()
    print("Extending CDSs according to coverage...")
    extended_df = extend_cds_regions(max_coverage_df, pd.concat([positive_bedgraph_df, negative_bedgraph_df]))
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    extended_df.to_csv(outPath, sep = "\t", header=False, index=False)
    #print(max_coverage_df)
