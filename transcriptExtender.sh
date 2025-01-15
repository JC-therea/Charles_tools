#!/bin/bash

# Verifica si los argumentos necesarios fueron proporcionados
if [ "$#" -ne 2 ]; then
    echo "Use: $0 <input.gff> <out directory>"
    exit 1
fi

# Variables de entrada
input_gff="$1"
output_dir="$2"

# Verifica si el archivo de entrada existe
if [ ! -f "$input_gff" ]; then
    echo "Error: File '$input_gff' doesn't exist."
    exit 1
fi

# Crea el directorio de salida si no existe
mkdir -p "$output_dir"

# Nombre del archivo de salida
output_gff="${output_dir}/$(basename "$input_gff" .gff)_extended.gff"
output_gtf_agat="${output_dir}/$(basename "$input_gff" .gff)_extended.agat.gtf"

# Procesar el archivo
awk 'BEGIN {OFS="\t"}
    $3 == "gene" || $3 == "transcript" || $3 == "mRNA" {
        $4 = ($4 - 50 < 1) ? 1 : $4 - 50;  # Extiende el start 50 nucleótidos, asegurándote de que no sea menor a 1
        $5 = $5 + 50;                      # Extiende el end 50 nucleótidos
    }
    { print }' "$input_gff" > "$output_gff"

echo "Fixing exon ends and defining UTR regions with agat..."

# Now fix with agat

source ~/mambaforge/bin/activate agatTest

agat_convert_sp_gff2gtf.pl --gff $output_gff -o $output_gtf_agat > textFromAgat.out
echo "output file: $output_gtf_agat"
