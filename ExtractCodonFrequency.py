import sys
from Bio import SeqIO
import pandas as pd

try:
    functionName, file, out = sys.argv
except:
    print("error: ORFSearcher.py <Input nucleotide fasta file> <Output prefix>")
    quit()


#file = "/home/jmontanes/Documents/0-Important_files/Insecta/Outputs/DrosophilaMelanogasterCDS/Dmel_CDSs_N1.fa"
#out = "test.tsv"

codonDict = {"AAA":0, "AAT":0, "AAC":0, "AAG":0,
"ATA":0, "ATT":0, "ATC":0, "ATG":0,
"ACA":0, "ACT":0, "ACC":0, "ACG":0,
"AGA":0, "AGT":0, "AGC":0, "AGG":0,
"CAA":0, "CAT":0, "CAC":0, "CAG":0,
"CTA":0, "CTT":0, "CTC":0, "CTG":0,
"CCA":0, "CCT":0, "CCC":0, "CCG":0,
"CGA":0, "CGT":0, "CGC":0, "CGG":0,
"TAA":0, "TAT":0, "TAC":0, "TAG":0,
"TTA":0, "TTT":0, "TTC":0, "TTG":0,
"TCA":0, "TCT":0, "TCC":0, "TCG":0,
"TGA":0, "TGT":0, "TGC":0, "TGG":0,
"GAA":0, "GAT":0, "GAC":0, "GAG":0,
"GTA":0, "GTT":0, "GTC":0, "GTG":0,
"GCA":0, "GCT":0, "GCC":0, "GCG":0,
"GGA":0, "GGT":0, "GGC":0, "GGG":0
             }

def trinucleotideParser(nuclSeq):
    pos = 0
    seqLen = len(nuclSeq)
    codonDict = {"AAA": 0, "AAT": 0, "AAC": 0, "AAG": 0,
                 "ATA": 0, "ATT": 0, "ATC": 0, "ATG": 0,
                 "ACA": 0, "ACT": 0, "ACC": 0, "ACG": 0,
                 "AGA": 0, "AGT": 0, "AGC": 0, "AGG": 0,
                 "CAA": 0, "CAT": 0, "CAC": 0, "CAG": 0,
                 "CTA": 0, "CTT": 0, "CTC": 0, "CTG": 0,
                 "CCA": 0, "CCT": 0, "CCC": 0, "CCG": 0,
                 "CGA": 0, "CGT": 0, "CGC": 0, "CGG": 0,
                 "TAA": 0, "TAT": 0, "TAC": 0, "TAG": 0,
                 "TTA": 0, "TTT": 0, "TTC": 0, "TTG": 0,
                 "TCA": 0, "TCT": 0, "TCC": 0, "TCG": 0,
                 "TGA": 0, "TGT": 0, "TGC": 0, "TGG": 0,
                 "GAA": 0, "GAT": 0, "GAC": 0, "GAG": 0,
                 "GTA": 0, "GTT": 0, "GTC": 0, "GTG": 0,
                 "GCA": 0, "GCT": 0, "GCC": 0, "GCG": 0,
                 "GGA": 0, "GGT": 0, "GGC": 0, "GGG": 0
                 }
    while pos < seqLen:
        codon = nuclSeq[pos:(pos + 3)]
        codonDict[codon] += 1
        pos += 3
    return(codonDict)

# Get list of transcripts
transcriptID_list = []
sequences = SeqIO.index(file, "fasta")

for transcript in sequences:
    transcriptID_list.append(transcript)

dfpd = pd.DataFrame(index = codonDict.keys(), columns = transcriptID_list)

for transcript in sequences:
    codonDictTranscript = trinucleotideParser(str(sequences[transcript].seq).upper())
    dfpd[transcript] = dfpd.index.map(codonDictTranscript)

dfpd_noStops = dfpd.drop(["TAG","TGA","TAA"])

AbsolutFrequencies = dfpd_noStops.sum(axis=1)
RelativeFrequencies = AbsolutFrequencies / sum(AbsolutFrequencies)
dfpd_noStops["AbsolutFrequencies"] = AbsolutFrequencies
dfpd_noStops["RelativeFrequencies"] = RelativeFrequencies

dfpd_noStops.to_csv(out, sep="\t", index=True, header=True)
