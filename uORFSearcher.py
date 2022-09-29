%%writefile "Inputs/Scripts/uORFSearcher.py"

import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import glob

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

    return idxs

def find_orfs_in_frame(sequence, threshold, name):
    """ Finds all valid open reading frames in the string 'sequence' (in frame), and
    returns them as a list"""

    orf_dict = {}
    frames = [0,1,2]
    for frame in frames:

        starts_a = find_all_in_frame(sequence[frame:], 'ATG')
        starts_t = find_all_in_frame(sequence[frame:], 'TTG')
        starts_c = find_all_in_frame(sequence[frame:], 'CTG')
        starts_g = find_all_in_frame(sequence[frame:], 'GTG')
        stop_amber = find_all_in_frame(sequence[frame:], 'TAG')
        stop_ochre = find_all_in_frame(sequence[frame:], 'TAA')
        stop_umber = find_all_in_frame(sequence[frame:], 'TGA')
        starts = starts_a + starts_c + starts_g + starts_t
        stops = stop_amber + stop_ochre + stop_umber
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

                    orf_obj = orf_object(sequence[bigger_start:stop + 3], bigger_start, stop + 3, name)
                    orfs.append(orf_obj)
            else:
                possible_starts = [start for start in starts if start < stop and (start - stop) % 3 == 0 and start > last_stop]
                possible_starts_above_threshold = [start for start in possible_starts if
                                                   len(sequence[start:stop + 3]) >= int(threshold) * 3]
                last_stop = stop
                if len(possible_starts_above_threshold) > 0:
                    possible_starts_above_threshold.sort()
                    bigger_start = possible_starts_above_threshold[0]
                    orf_obj = orf_object(sequence[bigger_start:stop + 3], bigger_start, stop + 3, name)
                    orfs.append(orf_obj)

        orfs.sort(key=lambda x: x.start, reverse=False)
        orf_dict[frame] = orfs
    return orf_dict


def from_orfs_to_gff(orfDict_pos, orfDict_neg, fastaseq,firstORF = 0):
    colnames = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attr"]
    frames = [0, 1, 2]

    totalEntries_pos = sum([len(orfList) for orfList in orfDict_pos.values()])
    totalEntries_neg = sum([len(orfList) for orfList in orfDict_neg.values()])

    chr = orfDict_pos[0][0].seqName

    # + strand
    gffPos = pd.DataFrame("x", index = range(totalEntries_pos), columns=colnames)
    gffPos.chr = chr
    gffPos.source = "ORFSearcher"
    gffPos.feature = "ORF"
    gffPos.score = "."
    gffPos.frame = "."
    gffPos.strand = "+"

    startListOfListsPos = []
    endListOfListsPos = []
    for frame in frames:
        startListOfListsPos.append([o.start for o in orfDict_pos[frame]])
        endListOfListsPos.append([o.end for o in orfDict_pos[frame]])

    startListPos = sum(startListOfListsPos, [])
    endListPos = sum(endListOfListsPos, [])

    gffPos.start = startListPos
    gffPos.end = endListPos

    ini = firstORF
    lastORF = ini + len(startListPos)
    attr = ['ID=ORF-' + str(i).zfill(8) for i in range(ini, lastORF)]
    gffPos.attr = attr
    # - strand
    gffNeg = pd.DataFrame("x", index = range(totalEntries_neg), columns=colnames)
    gffNeg.chr = chr
    gffNeg.source = "ORFSearcher"
    gffNeg.feature = "ORF"
    gffNeg.score = "."
    gffNeg.frame = "."
    gffNeg.strand = "-"
    startListOfListsNeg = []
    endListOfListsNeg = []
    for frame in frames:
        endListOfListsNeg.append([o.start for o in orfDict_neg[frame]])
        startListOfListsNeg.append([o.end for o in orfDict_neg[frame]])

    startListNeg = len(fastaseq) - pd.array(sum(startListOfListsNeg, []))
    endListNeg = len(fastaseq) - pd.array(sum(endListOfListsNeg, []))
    gffNeg.start = startListNeg
    gffNeg.end = endListNeg
    ini = lastORF + 1
    lastORF = ini + len(startListNeg)
    attr = ['ID=ORF-' + str(i).zfill(8) for i in range(ini, lastORF)]
    gffNeg.attr = attr
    gff = pd.concat([gffPos, gffNeg])
    gff.index = range(len(startListPos) + len(startListNeg))

    return(gff)






######################################################
###### Read each file in the desired directory #######
######################################################

# Input
fasta = "Inputs/TestTRANSCRIPTS.fa"
# Read the input as fasta
fastaIO = SeqIO.index(fasta, "fasta")
output = "potato.gff"

# Here we change aa threshold
threshold = 16
#############################

for sequence in fastaIO:
    orfDict_Pos = find_orfs_in_frame(fastaIO[sequence].seq, threshold, sequence)
    orfDict_Neg = find_orfs_in_frame(fastaIO[sequence].seq, threshold, sequence)



for sequence in fastaIO:
    if sequence:
        orfs_1[sequence] = find_orfs_in_frame(fastaIO[sequence].seq, threshold)
for sequence in fastaIO:
    if sequence:
        orfs_2[sequence] = find_orfs_in_frame(fastaIO[sequence].seq[1:], threshold)
for sequence in fastaIO:
    if sequence:
        orfs_3[sequence] = find_orfs_in_frame(fastaIO[sequence].seq[2:], threshold)

print("Loading ORFs from " + yeast)
CDS_dictionary = {}

for transcript, ORF_objects in orfs_1.items():
    if transcript:
        longest_ORF_sequence = ''
        longest_ORF_length = 0
        n = 0
        for single_ORF_object in ORF_objects:
            n += 1
            longest_ORF_sequence = single_ORF_object.sequence
            transcript_add = transcript + "_" + str(n) + "-frame_1"
            CDS_dictionary[transcript_add] = longest_ORF_sequence
    else:
        pass

for transcript, ORF_objects in orfs_2.items():
    if transcript:
        longest_ORF_sequence = ''
        longest_ORF_length = 0
        n = 0
        for single_ORF_object in ORF_objects:
            n += 1
            longest_ORF_sequence = single_ORF_object.sequence
            transcript_add = transcript + "_" + str(n) + "-frame_2"
            CDS_dictionary[transcript_add] = longest_ORF_sequence
    else:
        pass

for transcript, ORF_objects in orfs_3.items():
    if transcript:
        longest_ORF_sequence = ''
        longest_ORF_length = 0
        n = 0
        for single_ORF_object in ORF_objects:
            n += 1
            longest_ORF_sequence = single_ORF_object.sequence
            transcript_add = transcript + "_" + str(n) + "-frame_3"
            CDS_dictionary[transcript_add] = longest_ORF_sequence
    else:
        pass
stored_ORFs = open(output, 'w+')
print("Wrote " + str(len(CDS_dictionary)) + " CDSs")

for sequence in CDS_dictionary:
    fasta_index = ">" + sequence + "\n"
    stored_ORFs.write(fasta_index)
    nucleotide_sequence = str(CDS_dictionary[sequence].translate()) + "\n"
    stored_ORFs.write(nucleotide_sequence)

stored_ORFs.close()


fastaFile = "Inputs/TestTRANSCRIPTS.fa"
fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
fasta = SeqIO.index(open(fastaFile),'fasta')

for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
