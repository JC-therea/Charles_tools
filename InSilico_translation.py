from Bio import SeqIO
import sys
from Bio.Seq import Seq
import glob


import sys

try:
    Path_to_files = sys.argv[1]
    AA_treshold = int(sys.argv[2])
except:
    print("error: inSilico_translation.py <Directory with multifasta files> <Amino acid threshold>")
    quit()

############# Define OBJECTS ##############
class orf_object:
    def __init__(self, sequence, start, end):
        self.sequence = sequence
        self.start = start
        self.end = end
        self.length = len(sequence)

###########################################

############# Define FUNCTIONS ############

def find_all_in_frame(sequence, subsequence):
    ''' Returns a list of indexes within sequence that are the start of subsequence and are in frame'''
    start = 0
    idxs = []
    next_idx = sequence.find(subsequence, start)
    while next_idx != -1:
        if next_idx % 3 == 0:
            idxs.append(next_idx)
        start = next_idx + 1  # Move past this on the next time around
        next_idx = sequence.find(subsequence, start)

    return idxs

def find_orfs_in_frame(sequence, threshold):
    """ Finds all valid open reading frames in the string 'sequence' (in frame), and
    returns them as a list"""

    #First frame

    starts_a = find_all_in_frame(sequence, 'ATG')
    starts_t = find_all_in_frame(sequence, 'TTG')
    starts_c = find_all_in_frame(sequence, 'CTG')
    starts_g = find_all_in_frame(sequence, 'GTG')
    stop_amber = find_all_in_frame(sequence, 'TAG')
    stop_ochre = find_all_in_frame(sequence, 'TAA')
    stop_umber = find_all_in_frame(sequence, 'TGA')
    starts = starts_a #+ starts_c + starts_g + starts_t
    stops = stop_amber + stop_ochre + stop_umber
    starts.sort()
    stops.sort()

    orfs = []

    last_stop = -1
    for stop in stops:
        if stop == stops[0]:
            last_stop = -1
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

                orf_obj = orf_object(sequence[bigger_start:stop + 3], bigger_start, stop + 3)
                orfs.append(orf_obj)
        else:
            possible_starts = [start for start in starts if start < stop and (start - stop) % 3 == 0 and start > last_stop]
            possible_starts_above_threshold = [start for start in possible_starts if
                                               len(sequence[start:stop + 3]) >= int(threshold) * 3]
            last_stop = stop
            if len(possible_starts_above_threshold) > 0:
                possible_starts_above_threshold.sort()
                bigger_start = possible_starts_above_threshold[0]
                orf_obj = orf_object(sequence[bigger_start:stop + 3], bigger_start, stop + 3)
                orfs.append(orf_obj)

    orfs.sort(key=lambda x: len(x.sequence), reverse=True)
    return orfs  # returns the locations of all ORFs longer than threshold in the fasta sequence

###########################################################################################

###### Read each file in the desired directory #######

files = [f for f in glob.glob(Path_to_files + "/*.fa", recursive=True)]

for f in files:

    # Input
    fasta = f
    # Read the input as fasta
    cdna = SeqIO.index(fasta, "fasta")
    # Build the output path
    output = f.replace(".fa",'_AA_predicted.fa')

    # One ORF per frame
    orfs_1 = {}
    orfs_2 = {}
    orfs_3 = {}

    # Here we change aa threshold
    threshold = AA_treshold
    #############################

    for sequence in cdna:
        if sequence:
            orfs_1[sequence] = find_orfs_in_frame(cdna[sequence].seq, threshold)
    for sequence in cdna:
        if sequence:
            orfs_2[sequence] = find_orfs_in_frame(cdna[sequence].seq[1:], threshold)
    for sequence in cdna:
        if sequence:
            orfs_3[sequence] = find_orfs_in_frame(cdna[sequence].seq[2:], threshold)

    print("Loading ORFs from " + f)
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
    print(str(len(CDS_dictionary)) + " peptides have been written")


    for sequence in CDS_dictionary:
        fasta_index = ">"+sequence+"\n"
        stored_ORFs.write(fasta_index)
        nucleotide_sequence = str(CDS_dictionary[sequence].translate())+"\n"
        stored_ORFs.write(nucleotide_sequence)

    stored_ORFs.close()
