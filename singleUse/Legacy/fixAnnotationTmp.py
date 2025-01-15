import re

def filter_attributes(attribute):
    # Definir los patrones que queremos mantener
    patterns = ["ID=", "Parent=", "Name=", "Alias="]
    
    # Separar el atributo por ";"
    fields = attribute.split(";")
    
    # Filtrar los campos que comienzan con los patrones definidos
    filtered_fields = [field for field in fields if any(field.startswith(p) for p in patterns)]
    
    # Volver a unir los campos filtrados
    return ";".join(filtered_fields)

def process_gff(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Saltar las líneas de comentarios
            if line.startswith("#"):
                outfile.write(line)
                continue
            
            # Dividir la línea en columnas
            columns = line.strip().split("\t")
            
            # Filtrar los atributos (última columna)
            columns[-1] = filter_attributes(columns[-1])
            
            # Escribir la línea modificada en el archivo de salida
            outfile.write("\t".join(columns) + "\n")

# Llamada a la función para procesar el archivo
input_file = "/home/jmontanes/Documents/EvolutionNanopore/tmp3"
output_file = "/datasets/Common/ReferenceGenomes/Yeasts/CER/R64/Scer_R64.jmontanes.adapted.gff"
process_gff(input_file, output_file)
