from Bio import SeqIO

import sys

try:
    functionName, inputFile, outDir, temperature, pH, ion = sys.argv
except:
    print("error: tangoFormatter.py <inputFile.fa> <outDir> <Temperature (Kelvins)> <pH> <ionic strength>" )
    quit()

input_file = inputFile
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
i = 0
file = 1
outFile = inputFile.split("/")[-1].split(".fa")[0]
output_path = outDir  + outFile + "_TANGO_" + str(file) + ".tsv"
output = open(output_path,"w+")

for fasta in fasta_sequences:
    i+= 1
    if i == 950:
        output.close()
        file += 1
        output_path = outDir  + outFile + "_TANGO_" + str(file) + ".tsv"
        output = open(output_path,"w+")
        i = 0
    
    name, sequence = fasta.id, str(fasta.seq)
    line = name +" "+ "N" +" "+ "N" +" "+ pH +" "+ temperature +" "+ ion +" "+ sequence.split(".")[0] + "\n"
    output.writelines(line)

output.close()
