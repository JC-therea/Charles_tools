
# Run on pythonPackages environment

import sys
import pandas as pd
import warnings
warnings.filterwarnings("ignore")


try:
    functionName, file, output = sys.argv
except:
    print("error: AnnotationQC.py <Input file> <Output file>")
    #quit()

#file = "/home/jmontanes/Documents/EvolutionNanopore/Inputs/Park_2014/editing/Scer_Park2014.gtf"
#output = "Result.tsv"

colnames = ["chr","source","feature","start","end","score","strand","frame","attr"]
transcriptTypes = ["mRNA","ncRNA","transcript","snoRNA","snRNA","rRNA","pseudogene","antisense_RNA","tRNA","RNase_P_RNA","RNase_MRP_RNA","SRP_RNA"]
UTRTypes = ["five_prime_UTR","three_prime_UTR","five_prime_utr","three_prime_utr"]

if file.split(".")[-1] == "gff" or file.split(".")[-1] == "gff3":
    
    IDreplacePre = "ID="
    IDreplacePost = ";.*"
    genereplacePre = ".*Parent="
    genereplacePost = ";.*"
    
elif file.split(".")[-1] == "gtf":

    IDreplacePre = '.*transcript_id "'
    IDreplacePost = '";.*'
    genereplacePre = 'gene_id "'
    genereplacePost = '";.*'

else:
    
    print("error: Non valid file")
    quit()
    
    
df = pd.read_csv(file, sep = "\t", comment = "#", names = colnames, header = None)
transcripts = df[df.feature.isin(transcriptTypes)]
transcripts['ID'] = transcripts["attr"].str.replace(IDreplacePre,'').str.replace(IDreplacePost,"")
transcripts['gene'] = transcripts["attr"].str.replace(genereplacePre,'').str.replace(genereplacePost,"")
mRNASize = transcripts.end.subtract(transcripts.start) + 1 
mRNASize = mRNASize.set_axis(transcripts.ID)
mRNASize.name = "mRNA_size"
# Calculate CDS region
CDS = df[df.feature.isin(["CDS"])]
CDS['ID'] = CDS["attr"].str.replace(IDreplacePre,'').str.replace(IDreplacePost,"")
CDS['gene'] = CDS["attr"].str.replace(genereplacePre,'').str.replace(genereplacePost,"")

CDS["Size"] = CDS.end.subtract(CDS.start) + 1 
CDSsum = CDS.groupby('gene')['Size'].transform("sum")
CDSSize = CDSsum.set_axis(CDS.ID)
CDSSize.name = "CDS_size"
CDSSize = CDSSize[~CDSSize.index.duplicated(keep='first')]
#########################################################################
# Calculate exon size
exon = df[df.feature.isin(["exon"])]
exon['ID'] = exon["attr"].str.replace(IDreplacePre,'').str.replace(IDreplacePost,"")
exon['gene'] = exon["attr"].str.replace(genereplacePre,'').str.replace(genereplacePost,"")

exon["Size"] = exon.end.subtract(exon.start) + 1 
if file.split(".")[-1] == "gff" or file.split(".")[-1] == "gff3":
    
    exonsum = exon.groupby('gene')['Size'].transform("sum")
    exonSize = exonsum.set_axis(exon.ID)
    exonSize.name = "Exon_size"
    exonSize = exonSize[~exonSize.index.duplicated(keep='first')]
    
else:
    
    exonsum = exon.groupby('ID')['Size'].transform("sum")
    exonSize = exonsum.set_axis(exon.ID)
    exonSize.name = "Exon_size"
    exonSize = exonSize[~exonSize.index.duplicated(keep='first')]

    
# Check if there is some UTR annotated or not
if len(df[df.feature.isin(UTRTypes)]) == 0:
    UTRSeries = pd.Series([0] * len(transcripts), name = "UTR_size")
    UTRSeries = UTRSeries.set_axis(transcripts.ID)
else:
    UTRs = df[df.feature.isin(UTRTypes)]
    UTRs['ID'] = UTRs["attr"].str.replace(IDreplacePre,'').str.replace(IDreplacePost,"")
    UTRs['gene'] = UTRs["attr"].str.replace(genereplacePre,'').str.replace(genereplacePost,"")
    UTRs["Size"] = UTRs.end.subtract(UTRs.start) + 1
    UTRsum = UTRs.groupby('gene')['Size'].transform("sum")
    UTRSeries = UTRsum.set_axis(UTRs.ID)
    UTRSeries.name = "UTR_size"
    UTRSeries = UTRSeries[~UTRSeries.index.duplicated(keep='first')]
# Split by 5' and 3' UTR
if len(df[df.feature.isin(["five_prime_UTR", "five_prime_utr"])]) == 0:
    UTR5Series = pd.Series([0] * len(transcripts), name = "UTR5_size")
    UTR5Series = UTR5Series.set_axis(transcripts.ID)
else:
    UTR5s = df[df.feature.isin(["five_prime_UTR", "five_prime_utr"])]
    UTR5s['ID'] = UTR5s["attr"].str.replace(IDreplacePre,'').str.replace(IDreplacePost,"")
    UTR5s['gene'] = UTR5s["attr"].str.replace(genereplacePre,'').str.replace(genereplacePost,"")
    UTR5s["Size"] = UTR5s.end.subtract(UTR5s.start) + 1
    UTR5sum = UTR5s.groupby('gene')['Size'].transform("sum")
    UTR5Series = UTR5sum.set_axis(UTR5s.ID)
    UTR5Series.name = "UTR5_size"
    UTR5Series = UTR5Series[~UTR5Series.index.duplicated(keep='first')]

if len(df[df.feature.isin(["three_prime_UTR", "three_prime_utr"])]) == 0:
    UTR3Series = pd.Series([0] * len(transcripts), name = "UTR3_size")
    UTR3Series = UTR3Series.set_axis(transcripts.ID)
else:
    UTR3s = df[df.feature.isin(["three_prime_UTR", "three_prime_utr"])]
    UTR3s['ID'] = UTR3s["attr"].str.replace(IDreplacePre,'').str.replace(IDreplacePost,"")
    UTR3s['gene'] = UTR3s["attr"].str.replace(genereplacePre,'').str.replace(genereplacePost,"")
    UTR3s["Size"] = UTR3s.end.subtract(UTR3s.start) + 1
    UTR3sum = UTR3s.groupby('gene')['Size'].transform("sum")
    UTR3Series = UTR3sum.set_axis(UTR3s.ID)
    UTR3Series.name = "UTR3_size"
    UTR3Series = UTR3Series[~UTR3Series.index.duplicated(keep='first')]

#########################################################################
frame = { 'Transcripts': transcripts['ID'],"Transcript_type": transcripts['feature'], 'Genes': transcripts['gene']}

exportFile = pd.DataFrame(data = frame)
exportFile = exportFile.set_index('Transcripts').copy()
exportFile = exportFile.merge(mRNASize,left_index=True, right_index=True).copy()
exportFile = exportFile.merge(exonSize,left_index=True, right_index=True, how = "outer").copy()
exportFile = exportFile.merge(CDSSize,left_index=True, right_index=True, how = "outer").copy()
exportFile = exportFile.merge(UTRSeries,left_index=True, right_index=True, how = "outer").copy()
exportFile = exportFile.merge(UTR5Series,left_index=True, right_index=True, how = "outer").copy()
exportFile = exportFile.merge(UTR3Series,left_index=True, right_index=True, how = "outer").copy()
exportFile = exportFile.fillna(0).copy()
exportFile["Intron_size"] = exportFile.mRNA_size.subtract(exportFile.Exon_size)
exportFile.index.name = "Transcript"
exportFile.to_csv(output, sep = "\t")
    


