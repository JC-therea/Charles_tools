import sys
import pandas as pd

try:
    functionName, VCF_path, output = sys.argv
except:
    print("error: Calculate_transition_matrix.py <VCF file path> <Output tsv file>")
    #quit()

# VCF_path = "Parsing_SNPs/SNP_files/Light_SN_SNP_AF_0_05_Scer.gvcf"
Transition_matrix_full = output + "_full_info.tsv"
Transition_matrix = output + ".tsv"

SNP_table = pd.read_table(VCF_path, sep='\t', comment="#", dtype='str', header=None)
print("Transform column 2 into numeric...")
SNP_table_v2 = SNP_table.rename(columns={0: "Chr", 1: "Position", 3:"Reference", 4:"Alternative"})
SNP_table_v2["Position"] = pd.to_numeric(SNP_table_v2["Position"])

AA = 0
AC = sum((SNP_table_v2["Reference"] == "A") & (SNP_table_v2["Alternative"] == "C"))
AG = sum((SNP_table_v2["Reference"] == "A") & (SNP_table_v2["Alternative"] == "G"))
AT = sum((SNP_table_v2["Reference"] == "A") & (SNP_table_v2["Alternative"] == "T"))

CA = sum((SNP_table_v2["Reference"] == "C") & (SNP_table_v2["Alternative"] == "A"))
CC = 0
CG = sum((SNP_table_v2["Reference"] == "C") & (SNP_table_v2["Alternative"] == "G"))
CT = sum((SNP_table_v2["Reference"] == "C") & (SNP_table_v2["Alternative"] == "T"))

GA = sum((SNP_table_v2["Reference"] == "G") & (SNP_table_v2["Alternative"] == "A"))
GC = sum((SNP_table_v2["Reference"] == "G") & (SNP_table_v2["Alternative"] == "C"))
GG = 0
GT = sum((SNP_table_v2["Reference"] == "G") & (SNP_table_v2["Alternative"] == "T"))

TA = sum((SNP_table_v2["Reference"] == "T") & (SNP_table_v2["Alternative"] == "A"))
TC = sum((SNP_table_v2["Reference"] == "T") & (SNP_table_v2["Alternative"] == "C"))
TG = sum((SNP_table_v2["Reference"] == "T") & (SNP_table_v2["Alternative"] == "G"))
TT = 0

expected = pd.DataFrame({'A': [AA,CA,GA,TA],
                         'C': [AC,CC,GC,TC],
                         'G': [AG,CG,GG,TG],
                         'T': [AT,CT,GT,TT]},
                         index = ['A','C','G','T'])


print(expected)

print("A -> C")
print(AC)

print("T -> G")
print(TG)

expected.to_csv(Transition_matrix_full, sep = "\t", index = ['A','C','G','T'], header = ['A','C','G','T'])
expected.to_csv(Transition_matrix, sep = "\t", header = False, index = False)