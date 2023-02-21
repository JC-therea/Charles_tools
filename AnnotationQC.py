
# Run on pythonPackages environment

import sys
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# It has to be fixed por gff files

try:
    functionName, file, output = sys.argv
except:
    print("error: AnnotationQC.py <Input file> <Output file>")
    #quit()

#file = "/projects_eg/projects/jmontanes/EvolutionNanopore/Annotation/Tests/UpdateOnly-WorkWithChapuza/ScerBloom/update_results/Saccharomyces_cerevisiaeJC.gff3"
#output = "Result.tsv"

colnames = ["chr","source","feature","start","end","score","strand","frame","attr"]
transcriptTypes = ["mRNA","ncRNA","transcript","snoRNA","snRNA","rRNA","pseudogene","antisense_RNA","tRNA","RNase_P_RNA","RNase_MRP_RNA","SRP_RNA"]
UTRTypes = ["five_prime_UTR","three_prime_UTR","five_prime_utr","three_prime_utr"]

if file.split(".")[-1] == "gff" or file.split(".")[-1] == "gff3":
    
    IDreplacePre = "ID="
    IDreplacePost = ";.*"
    genereplacePre = ".*Parent="
    genereplacePost = ";.*"

    IDreplacePreNoTransc = ".*Parent="
    IDreplacePostNoTransc = ";.*"
    genereplacePreNoTransc = ".*Parent="
    genereplacePostNoTransc = "-T.;.*"
    
elif file.split(".")[-1] == "gtf":

    IDreplacePre = '.*transcript_id "'
    IDreplacePost = '";.*'
    genereplacePre = 'gene_id "'
    genereplacePost = '";.*'

    IDreplacePreNoTransc = '.*transcript_id "'
    IDreplacePostNoTransc = '";.*'
    genereplacePreNoTransc = 'gene_id "'
    genereplacePostNoTransc = '";.*'

else:
    
    print("error: Non valid file")
    quit()
    
    
df = pd.read_csv(file, sep = "\t", comment = "#", names = colnames, header = None)
transcripts = df[df.feature.isin(transcriptTypes)]
transcripts = transcripts.assign(ID = transcripts["attr"].str.replace(IDreplacePre,'',regex=True).str.replace(IDreplacePost,"",regex=True))
transcripts = transcripts.assign(gene = transcripts["attr"].str.replace(genereplacePre,'',regex = True).str.replace(genereplacePost,"", regex = True))

mRNASize = transcripts.end.subtract(transcripts.start)
mRNASize = mRNASize.set_axis(transcripts.ID)
mRNASize.name = "mRNA_size"
# Calculate CDS region
CDS = df[df.feature.isin(["CDS"])]
CDS = CDS.assign(ID=CDS["attr"].str.replace(IDreplacePreNoTransc, '', regex=True).str.replace(IDreplacePostNoTransc, "", regex=True))
CDS = CDS.assign(gene=CDS["attr"].str.replace(genereplacePreNoTransc, '', regex=True).str.replace(genereplacePostNoTransc, "",regex=True))
CDS["Size"] = CDS.end.subtract(CDS.start)
CDSsum = CDS.groupby('ID')['Size'].transform("sum")
CDSSize = CDSsum.set_axis(CDS.ID)
CDSSize.name = "CDS_size"
CDSSize = CDSSize[~CDSSize.index.duplicated(keep='first')]

#########################################################################
# Calculate exon size
exon = df[df.feature.isin(["exon"])]
exon = exon.assign(ID=exon["attr"].str.replace(IDreplacePreNoTransc, '', regex=True).str.replace(IDreplacePostNoTransc, "", regex=True))
exon = exon.assign(gene=exon["attr"].str.replace(genereplacePreNoTransc, '', regex=True).str.replace(genereplacePostNoTransc, "",regex=True))
exon["Size"] = exon.end.subtract(exon.start)

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
    UTRs = UTRs.assign(ID=UTRs["attr"].str.replace(IDreplacePreNoTransc, '', regex=True).str.replace(IDreplacePostNoTransc, "", regex=True))
    UTRs = UTRs.assign(gene=UTRs["attr"].str.replace(genereplacePreNoTransc, '', regex=True).str.replace(genereplacePostNoTransc, "",regex=True))
    UTRs["Size"] = UTRs.end.subtract(UTRs.start)
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
    UTR5s = UTR5s.assign(ID=UTRs["attr"].str.replace(IDreplacePreNoTransc, '', regex=True).str.replace(IDreplacePostNoTransc, "", regex=True))
    UTR5s = UTR5s.assign(gene=UTRs["attr"].str.replace(genereplacePreNoTransc, '', regex=True).str.replace(genereplacePostNoTransc, "",regex=True))
    UTR5s["Size"] = UTR5s.end.subtract(UTR5s.start)
    UTR5sum = UTR5s.groupby('gene')['Size'].transform("sum")
    UTR5Series = UTR5sum.set_axis(UTR5s.ID)
    UTR5Series.name = "UTR5_size"
    UTR5Series = UTR5Series[~UTR5Series.index.duplicated(keep='first')]

if len(df[df.feature.isin(["three_prime_UTR", "three_prime_utr"])]) == 0:
    UTR3Series = pd.Series([0] * len(transcripts), name = "UTR3_size")
    UTR3Series = UTR3Series.set_axis(transcripts.ID)
else:
    UTR3s = df[df.feature.isin(["three_prime_UTR", "three_prime_utr"])]
    UTR3s = UTR3s.assign(ID=UTRs["attr"].str.replace(IDreplacePreNoTransc, '', regex=True).str.replace(IDreplacePostNoTransc, "", regex=True))
    UTR3s = UTR3s.assign(gene=UTRs["attr"].str.replace(genereplacePreNoTransc, '', regex=True).str.replace(genereplacePostNoTransc, "",regex=True))
    UTR3s["Size"] = UTR3s.end.subtract(UTR3s.start)
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
    


