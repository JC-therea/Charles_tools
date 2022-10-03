
import sys
from Bio import SeqIO
import pandas as pd

try:
    functionName, file, AA_min, mode, outputPrefix = sys.argv
except:
    print("error: ORFSearcher.py <Input file> <Minimum ORF size in AA> <Mode, available options 'All' 'ATG' and 'Stop'> <Output prefix>")
    quit()

modeOptions = ['All', 'ATG', 'Stop']

if mode not in mode:
    print("Mode is not well defined mode set as All")
    mode = "All"

############# Define OBJECTS ##############
class orf_object:
    def __init__(self, sequence, start, end, seqName):
        self.sequence = sequence
        self.start = start
        self.end = end
        self.length = len(sequence)
        self.seqName = str(seqName)

###########################################

############# Define FUNCTIONS ############

def find_all_in_frame(sequence, subsequence):
    #''' Returns a list of indexes within sequence that are the start of subsequence and are in frame'''
    start = 0
    idxs = []
    next_idx = sequence.find(subsequence, start)
    while next_idx != -1:
        if next_idx % 3 == 0:
            idxs.append(next_idx)
        start = next_idx + 1  # Move past this on the next time around
        next_idx = sequence.find(subsequence, start)

    idxs_1_based = [x+1 for x in idxs]
    return idxs_1_based

def find_orfs_in_frame(sequence, threshold, name, mode):
    """ Finds all valid open reading frames in the string 'sequence' (in frame), and
    returns them as a list"""

    orf_dict = {}
    frames = [0,1,2]
    for frame in frames:

        if mode == "All":
            starts_a = find_all_in_frame(sequence[frame:], 'ATG')
            starts_t = find_all_in_frame(sequence[frame:], 'TTG')
            starts_c = find_all_in_frame(sequence[frame:], 'CTG')
            starts_g = find_all_in_frame(sequence[frame:], 'GTG')

            starts_a_gff_adjustment = [x + frame for x in starts_a]
            starts_t_gff_adjustment = [x + frame for x in starts_t]
            starts_c_gff_adjustment = [x + frame for x in starts_c]
            starts_g_gff_adjustment = [x + frame for x in starts_g]

            starts = starts_a_gff_adjustment + starts_t_gff_adjustment + starts_c_gff_adjustment + starts_g_gff_adjustment
        elif mode == "ATG":
            starts_a = find_all_in_frame(sequence[frame:], 'ATG')
            starts_a_gff_adjustment = [x + frame for x in starts_a]
            starts = starts_a_gff_adjustment
            
        
		# Bug with mode Stop
		# This part seems to work 
        elif mode == "Stop":
            stop_amber = find_all_in_frame(sequence[frame:], 'TAG')
            stop_ochre = find_all_in_frame(sequence[frame:], 'TAA')
            stop_umber = find_all_in_frame(sequence[frame:], 'TGA')
            stop_amber_gff_adjustment = [x + frame for x in stop_amber]
            stop_ochre_gff_adjustment = [x + frame for x in stop_ochre]
            stop_umber_gff_adjustment = [x + frame for x in stop_umber]
            starts = stop_amber_gff_adjustment + stop_ochre_gff_adjustment + stop_umber_gff_adjustment

        stop_amber = find_all_in_frame(sequence[frame:], 'TAG')
        stop_ochre = find_all_in_frame(sequence[frame:], 'TAA')
        stop_umber = find_all_in_frame(sequence[frame:], 'TGA')

        stop_amber_gff_adjustment = [x + frame for x in stop_amber]
        stop_ochre_gff_adjustment = [x + frame for x in stop_ochre]
        stop_umber_gff_adjustment = [x + frame for x in stop_umber]

        stops = stop_amber_gff_adjustment + stop_ochre_gff_adjustment + stop_umber_gff_adjustment
        starts.sort()
        stops.sort()

        orfs = []

        last_stop = 0
        for stop in stops:
            if stop == stops[0]:
                last_stop = 0
                # Start is before the stop codon and after the last stop codon
                possible_starts = [start for start in starts if start < stop and (start - stop) % 3 == 0 and start > last_stop]
                # Limit size by the threshold
                possible_starts_above_threshold = [start for start in possible_starts if len(sequence[start:stop + 3]) >= int(threshold) * 3]
                # Update last stop codon
                last_stop = stop
                # If there is some start that in possible_starts_above_threshold then:
                if len(possible_starts_above_threshold) > 0:
                    # Sort the start positions
                    possible_starts_above_threshold.sort()
                    # Choose the bigger CDS
                    bigger_start = possible_starts_above_threshold[0]

                    orf_obj = orf_object(sequence[bigger_start:stop + 2], bigger_start, stop + 2, name)
                    orfs.append(orf_obj)
            else:
                possible_starts = [start for start in starts if start < stop and (start - stop) % 3 == 0 and start > last_stop]
                possible_starts_above_threshold = [start for start in possible_starts if
                                                   len(sequence[start:stop + 3]) >= int(threshold) * 3]
                last_stop = stop
                if len(possible_starts_above_threshold) > 0:
                    possible_starts_above_threshold.sort()
                    bigger_start = possible_starts_above_threshold[0]
                    orf_obj = orf_object(sequence[bigger_start:stop + 2], bigger_start, stop + 2, name)
                    orfs.append(orf_obj)

        orfs.sort(key=lambda x: x.start, reverse=False)
        orf_dict[frame] = orfs
    return orf_dict

# This part does not work with mode == Stop
def from_orfs_to_gff(orfDict_pos, orfDict_neg, fastaseq,firstORF = 0):

    frames = [0, 1, 2]
    chr = orfDict_pos[0][0].seqName

    Start = []
    End = []
    AttrList = []

    # + strand
    startListOfListsPos = []
    endListOfListsPos = []
    for frame in frames:
        startListOfListsPos.append([o.start for o in orfDict_pos[frame]])
        endListOfListsPos.append([o.end for o in orfDict_pos[frame]])

    startListPos = sum(startListOfListsPos, [])
    endListPos = sum(endListOfListsPos, [])

    Start.extend(startListPos)
    End.extend(endListPos)

    ini = firstORF
    lastORF = ini + len(startListPos)
    attr = ['ID=ORF-' + str(i).zfill(8) for i in range(ini, lastORF)]
    AttrList.extend(attr)

    # - strand
    startListOfListsNeg = []
    endListOfListsNeg = []
    for frame in frames:
        endListOfListsNeg.append([o.start for o in orfDict_neg[frame]])
        startListOfListsNeg.append([o.end for o in orfDict_neg[frame]])

    startListNeg = len(fastaseq) - pd.array(sum(startListOfListsNeg, [])) + 1
    endListNeg = len(fastaseq) - pd.array(sum(endListOfListsNeg, [])) + 1

    Start.extend(startListNeg)
    End.extend(endListNeg)

    ini = lastORF
    lastORF = ini + len(startListNeg)
    attr = ['ID=ORF-' + str(i).zfill(8) for i in range(ini, lastORF)]
    AttrList.extend(attr)

    chrList = [chr] * len(Start)
    source = ["ORFSearcher"] * len(Start)
    feature = ["ORF"] * len(Start)
    score = ["."] * len(Start)
    frame = ["."] * len(Start)
    Strand = ["+"] * len(startListPos) + ["-"] * len(startListNeg)

    return(chrList, source, feature, Start, End, score, Strand, frame, AttrList, lastORF)


######################################################
###### Read each file in the desired directory #######
######################################################

# Input
fasta = file
outFile = outputPrefix + "_ORFs.gff3"
# Read the input as fasta
fastaIO = SeqIO.index(fasta, "fasta")
output = "potato.gff"

# Here we change aa threshold
threshold = AA_min
#############################
colnames = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attr"]

chrList = []
sourceList = []
featureList = []
startList = []
endList = []
scoreList = []
strandList = []
frameList = []
attrList = []
lastORF = 0
for sequence in fastaIO:
    print("Starting ORF search in " + sequence + "...", flush=True)
    orfDict_Pos = find_orfs_in_frame(str(fastaIO[sequence].seq).upper(), threshold, sequence, mode)
    print("50 % done in " + sequence, flush=True)
    orfDict_Neg = find_orfs_in_frame(str(fastaIO[sequence].seq.reverse_complement()).upper(), threshold, sequence, mode)
    print("All ORFs in " + sequence + " found, finishing last details...", flush=True)
    chr, source, feature, start, end, score, strand, frame, attr, lastORF = from_orfs_to_gff(orfDict_Pos,orfDict_Neg, fastaIO[sequence].seq, lastORF)
    chrList.extend(chr)
    sourceList.extend(source)
    featureList.extend(feature)
    startList.extend(start)
    endList.extend(end)
    scoreList.extend(score)
    strandList.extend(strand)
    frameList.extend(frame)
    attrList.extend(attr)
    print(sequence + " done!")

gff = pd.DataFrame("x", index = range(len(chrList)), columns=colnames)
gff.chr = chrList
gff.source = sourceList
gff.feature = featureList
gff.start = startList
gff.end = endList
gff.score = scoreList
gff.strand = strandList
gff.frame = frameList
gff.attr = attrList

gffSorted = gff.sort_values(["chr", "start"], ascending = [True, True])
gffSorted.attr = ['ID=ORF-' + str(i).zfill(8) for i in range(0, len(gff.attr))]
gffSorted.to_csv(outFile, sep="\t", index=False, header=False)
