import pandas as pd
import sys

# Cargar el archivo GTF
def load_gtf(file):
    gtf = pd.read_csv(file, sep='\t', comment='#', header=None,
                      names=["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    return gtf

# Función para extraer el transcript_id y gene_id
def extract_id(attribute, id_type):
    try:
        # Extrae el valor del id_type (transcript_id o gene_id)
        return [item.split('"')[1] for item in attribute.split("; ") if id_type in item][0]
    except IndexError:
        return None

# Añadir el feature transcript
def add_transcript_features(gtf):
    gtf['transcript_id'] = gtf['attribute'].apply(lambda x: extract_id(x, "transcript_id"))
    gtf['gene_id'] = gtf['attribute'].apply(lambda x: extract_id(x, "gene_id"))

    # Agrupar por transcript_id y calcular las posiciones mínimas y máximas
    transcript_group = gtf.groupby('transcript_id').agg({'start': 'min', 'end': 'max', 'chr': 'first', 'strand': 'first', 'gene_id': 'first'}).reset_index()

    # Crear una nueva tabla con el feature "transcript"
    transcripts = transcript_group.copy()
    transcripts["feature"] = "transcript"
    transcripts["source"] = "python"
    transcripts["score"] = "."
    transcripts["frame"] = "."
    transcripts["attribute"] = transcripts.apply(lambda row: f'gene_id "{row.gene_id}"; transcript_id "{row.transcript_id}"', axis=1)

    return transcripts[['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]

# Añadir el feature gene
def add_gene_features(gtf):
    # Agrupar por gene_id y calcular las posiciones mínimas y máximas
    gene_group = gtf.groupby('gene_id').agg({'start': 'min', 'end': 'max', 'chr': 'first', 'strand': 'first'}).reset_index()

    # Crear una nueva tabla con el feature "gene"
    genes = gene_group.copy()
    genes["feature"] = "gene"
    genes["source"] = "python"
    genes["score"] = "."
    genes["frame"] = "."
    genes["attribute"] = genes.apply(lambda row: f'gene_id "{row.gene_id}"', axis=1)

    return genes[['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]

# Guardar el nuevo archivo GTF
def save_gtf(gtf, output_file):
    gtf.to_csv(output_file, sep='\t', header=False, index=False)

# Main
if __name__ == "__main__":
    # Cargar archivo GTF
    
    if len(sys.argv) != 3:
        print("Usage: python adjust_gtf.py input.gtf output.gtf")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    gtf = load_gtf(input_file)
    
    # Añadir las líneas de transcript y gene
    transcripts = add_transcript_features(gtf)
    genes = add_gene_features(gtf)
    
    # Combinar todo el GTF
    full_gtf = pd.concat([gtf, transcripts, genes], ignore_index=True).sort_values(by=['chr', 'start'])

    # Guardar el archivo
    save_gtf(full_gtf, output_file)
