import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="This program calculates a neutral model of amino acid substitutions and it compares it with the observations",
								 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c", "--codons", default="0", type=str, help="path of the codon estimation file", required=True)
parser.add_argument("-ns", "--nucleotideSubstitution", default="0", type=str, help="nucleotide substitutions matrix in the species (typically from introns)", required=True)
parser.add_argument("-m", "--mode", default="noDir", type=str, help="Select the mode of the program, if the codons that you introduce are from triplets mode is 'Dir' if you are introducing only pairwise alignments the keep the default mode 'noDir'", required=True)
parser.add_argument("-i", "--inputFiles", nargs='+', help="Add different tsv files with the codon change frequency", required=False)
parser.add_argument("-o", "--outFile", default="", type=str, help="Name of the output file", required=False)

args = parser.parse_args()
file_cod = args.codons
#file_cod = "/home/jmontanes/Documents/IQtree_Gene_duplication/Insects/DropboxThings/gain_of_acidic_residues/gain_of_acidic_residues/neutral_changes/codons.txt"
file_align_list = args.inputFiles
mode = args.mode
#mode = "Dir"
#file_align1="/home/jmontanes/Documents/IQtree_Gene_duplication/Insects/DropboxThings/gain_of_acidic_residues/gain_of_acidic_residues/AlignmentChangesN1_Drosophila_sechelia.tsv"
#file_align2="/home/jmontanes/Documents/IQtree_Gene_duplication/Insects/DropboxThings/gain_of_acidic_residues/gain_of_acidic_residues/AlignmentChangesN1_Drosophila_simulans.tsv"
#file_align_list = [file_align1, file_align2]
aa={
"A" : 1, "C" : 2, "D" : 3, "E" : 4, "F" : 5, "G" : 6, "H" : 7, "I" : 8, "K" : 9, "L" : 10,
"M" : 11, "N" : 12, "P" : 13, "Q" : 14, "R" : 15, "S" : 16, "T" : 17, "V" : 18, "W" : 19, "Y" : 20
}

aa_type={
"A" : "nonPolar", "C" : "Polar",  "D" : "Acidic", "E" : "Acidic", "F" : "nonPolar", "G" : "nonPolar", "H" : "Basic", "I" : "nonPolar", "K" : "Basic", "L" : "nonPolar",
"M" : "nonPolar", "N" : "Polar", "P" : "nonPolar", "Q" : "Polar", "R" : "Basic", "S" : "Polar", "T" : "Polar", "V" : "nonPolar", "W" : "nonPolar", "Y" : "Polar"
}

genet={"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y",
       "TGT":"C", "TGC":"C","TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

#first argument: file with codon frequencies 
#second argument: file with aa substitution frequencies of pairwise alignments D.mel-D.sech
#third argument: file with aa substitution frequencies of pairwise alignments D.mel-D.sim

#for set of CONSERVED proteins: python3 expected_aa_changes_cons.py codons_CONS.txt AlignmentChangesConserved_Drosophila_sechelia.tsv AlignmentChangesConserved_Drosophila_simulans.tsv
#for set of N1 proteins: python3 expected_aa_changes.py codons.txt AlignmentChangesN1_Drosophila_sechelia.tsv AlignmentChangesN1_Drosophila_simulans.tsv

###### Get the dictionary from tsv file
transitionMatrix = pd.read_csv(args.nucleotideSubstitution, sep="\t", header=0, index_col=0)
#transitionMatrix = pd.read_csv("/home/jmontanes/Documents/0-Important_files/Insecta/PNPS/Transition_matrix_introns_full_info.tsv", sep="\t", header=0, index_col=0)
totalCases = transitionMatrix.sum().sum()

subs={
"aa" : 0,
"ac" : transitionMatrix["C"]["A"] / totalCases,
"ag" : transitionMatrix["G"]["A"] / totalCases,
"at" : transitionMatrix["T"]["A"] / totalCases,
"ca" : transitionMatrix["A"]["C"] / totalCases,
"cc" : 0,
"cg" : transitionMatrix["G"]["C"] / totalCases,
"ct": transitionMatrix["T"]["C"] / totalCases,
"ga" : transitionMatrix["A"]["G"] / totalCases,
"gc" : transitionMatrix["C"]["G"] / totalCases,
"gg" : 0,
"gt" : transitionMatrix["T"]["G"] / totalCases,
"ta" : transitionMatrix["A"]["T"] / totalCases,
"tc" : transitionMatrix["C"]["T"] / totalCases,
"tg" : transitionMatrix["G"]["T"] / totalCases,
"tt" : 0}

######


aa_list=('A','C','D','E','F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

# new hashes to be created
tri_freq={}
prob={}

# fill in hash tri_freq codon frequencies 

fh=open(file_cod,'r')  
text=fh.readlines()
for line in text:
	row = line.split('\n')[0].split("\t")
	tri_freq[row[0]] = row[1]

# initialize hash prob 
for i in range(len(aa_list)):
	for j in range(len(aa_list)):
		aa_pair=aa_list[i] + aa_list[j]
		prob[aa_pair]=0

			
#########################################################################################################
#### build neutral model of amino acid change frequencies                               
#########################################################################################################

# algorithm: we go through each codon and mutate it in all positions, then check if this is associated with a codon for another amino acid, if so we add the probability of change to the corresponding variable

for key,value in tri_freq.items():
	#### changes in first position of a codon ####
	if ((key[0]) == "A"):  ## A
		#A to C
		mut=subs["ac"]  ## A to C
		codon_to_searchfor="C"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#A to G
		mut=subs["ag"] ## A to G
		codon_to_searchfor="G"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#A to T
		mut=subs["at"] ## A to T
		codon_to_searchfor="T"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
					
	if ((key[0]) == "C"): ## C
		#C to A
		mut=subs["ca"] ## C to A
		codon_to_searchfor="A"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#C to G
		mut=subs["cg"] ## C to G
		codon_to_searchfor="G"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
					
		#C to T
		mut=subs["ct"] ## C to T
		codon_to_searchfor="T"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
					
	if ((key[0]) == "G"): ## G
		#G to A
		mut=subs["ga"] ## G to A
		codon_to_searchfor="A"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#G to C
		mut=subs["gc"] ## G to C
		codon_to_searchfor="C"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
					
		#G to T
		mut=subs["gt"] ## G to T
		codon_to_searchfor="T"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
								
	if ((key[0]) == "T"): ## T
		#T to A
		mut=subs["ta"] ## T to A
		codon_to_searchfor="A"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#T to C
		mut=subs["tc"] ## T to C
		codon_to_searchfor="C"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
					
		#T to G
		mut=subs["tg"] ## T to G
		codon_to_searchfor="G"+key[1]+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
	#### changes in second position of a codon ####
	if ((key[1]) == "A"):  ## A
		#A to C
		mut=subs["ac"]  ## A to C
		codon_to_searchfor=key[0]+"C"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#A to G
		mut=subs["ag"]  ## A to G
		codon_to_searchfor=key[0]+"G"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
								
		#A to T
		mut=subs["at"]  ## A to T
		codon_to_searchfor=key[0]+"T"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
							
	if ((key[1]) == "C"):  ## C
		#C to A
		mut=subs["ca"]  ## C to A
		codon_to_searchfor=key[0]+"A"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#C to G
		mut=subs["cg"]  ## C to G
		codon_to_searchfor=key[0]+"G"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
							
		#C to T
		mut=subs["ct"]  ## C to T
		codon_to_searchfor=key[0]+"T"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
						
	if ((key[1]) == "G"):  ## G
		#G to A
		mut=subs["ga"]  ## G to A
		codon_to_searchfor=key[0]+"A"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#G to C
		mut=subs["gc"]  ## G to C
		codon_to_searchfor=key[0]+"C"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
							
		#G to T
		mut=subs["gt"]  ## G to T
		codon_to_searchfor=key[0]+"T"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
											
	if ((key[1]) == "T"):  ## T
		#T to A
		mut=subs["ta"]  ## T to A
		codon_to_searchfor=key[0]+"A"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#T to C
		mut=subs["tc"]  ## T to C
		codon_to_searchfor=key[0]+"C"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
							
		#T to G
		mut=subs["tg"]  ## T to G
		codon_to_searchfor=key[0]+"G"+key[2]
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
					
	#### changes in third position of a codon ####
	if ((key[2]) == "A"):  ## A
		#A to C
		mut=subs["ac"]  ## A to C
		codon_to_searchfor=key[0]+key[1]+"C"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#A to G
		mut=subs["ag"]  ## A to G
		codon_to_searchfor=key[0]+key[1]+"G"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
			
		#A to T
		mut=subs["at"]  ## A to T
		codon_to_searchfor=key[0]+key[1]+"T"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
	if ((key[2]) == "C"):  ## C
		#C to A
		mut=subs["ca"]  ## C to A
		codon_to_searchfor=key[0]+key[1]+"A"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#C to G
		mut=subs["cg"]  ## C to G
		codon_to_searchfor=key[0]+key[1]+"G"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#C to T
		mut=subs["ct"]  ## C to T
		codon_to_searchfor=key[0]+key[1]+"T"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
	if ((key[2]) == "G"):  ## G
		#G to A
		mut=subs["ga"]  ## G to A
		codon_to_searchfor=key[0]+key[1]+"A"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#G to C
		mut=subs["gc"]  ## G to C
		codon_to_searchfor=key[0]+key[1]+"C"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#G to T
		mut=subs["gt"]  ## G to T
		codon_to_searchfor=key[0]+key[1]+"T"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
	if ((key[2]) == "T"):  ## T
		#T to A
		mut=subs["ta"]  ## T to A
		codon_to_searchfor=key[0]+key[1]+"A"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#T to C
		mut=subs["tc"]  ## T to C
		codon_to_searchfor=key[0]+key[1]+"C"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
				
		#T to G
		mut=subs["tg"]  ## T to G
		codon_to_searchfor=key[0]+key[1]+"G"
		for keys,values in genet.items():
			if(keys==codon_to_searchfor):
				pair=genet[key]+values 
				prob[pair] = prob[pair] + (float(value)*float(mut))
						
													
# make substitutions of the same amino acid to be 0
for AA in aa_list:
	pair = AA+AA
	prob[pair] = 0

# store the expected probabilities in a file (with direction of change)
outName = args.outFile + 'neutral_model_probabilities.tsv'
pd.DataFrame.from_dict(prob, orient="index").to_csv(outName, sep="\t", header=False)

if mode == "noDir":
	#########################################################################################################
	#### calculate expected probabilities (with no direction of change, AE=EA) and compare with observed ####
	#########################################################################################################
	probNoDir = prob.copy()
	sum_prob_neu=0
	for i in range(len(aa_list)):
		for j in range(len(aa_list)):
			aa_pair1=aa_list[i] + aa_list[j]
			aa_pair2=aa_list[j] + aa_list[i]
			#print(aa_pair1," ",aa_pair2)
			if (aa_pair1 in probNoDir.keys() and aa_pair2 in probNoDir.keys()):
				probNoDir[aa_pair1] = probNoDir[aa_pair1] + probNoDir[aa_pair2]# add probabilities
				sum_prob_neu = sum_prob_neu + probNoDir[aa_pair1]
				removed=probNoDir.pop(aa_pair2) # remove second pair from the hash

	# store the expected probabilities in a file (with no direction of change)
	outName = args.outFile + 'neutral_model_probabilities_nodirectionofchange.tsv'
	pd.DataFrame.from_dict(probNoDir, orient="index").to_csv(outName, sep="\t", header=False)

	# get the observed data from files (only changes in the neutral model, which are at most one nucleotide away)

	probDf = pd.DataFrame.from_dict(probNoDir, orient="index")
	scaling_factor_obs_dict = {}
	pandasDfDict = {}
	pandasDfNames = []

	for fileIn in file_align_list:

		#fileIn = file_align_list[0]
		fileName = fileIn.split("/")[-1].split(".")[0]
		pandasDfNames.append(fileName)
		fileInDf = pd.read_csv(fileIn, sep="\t", index_col=0).drop("Perf")
		fileInDf_prob = fileInDf[fileInDf.index.isin(probNoDir.keys())].copy()
		fileInDf_prob.loc[fileInDf_prob.index.isin(probDf.index[probDf[0] == 0]), "Hits"] = 0
		pandasDfDict[fileName] = fileInDf_prob
		# calculate scaling factors (we make the total number of changes equal in expected and observed)
		scaling_factor_obs = fileInDf_prob.Hits.sum() / sum_prob_neu
		scaling_factor_obs_dict[fileName] = scaling_factor_obs

	#fres=open(args.outFile,"w+")
	outFileDf = probDf.copy()#.pop(0)
	for fileName in pandasDfNames:

		outFileDf["expected_" + fileName] = outFileDf[0] * scaling_factor_obs_dict[fileName]
		outFileDf["observed_" + fileName] = pandasDfDict[fileName].Hits

	firstLetter = outFileDf.index.get_level_values(outFileDf.index.name).str[0].tolist()
	firstType = [aa_type[x] for x in firstLetter]

	secondLetter = outFileDf.index.get_level_values(outFileDf.index.name).str[1].tolist()
	secondType = [aa_type[x] for x in secondLetter]

	Type = []

	for s1, s2 in zip(firstType, secondType):
		if s1[0] < s2[0]:
			Type.append(s1 + "-" + s2)
		else:
			Type.append(s2 + "-" + s1)

	outFileDf["Type"] = Type
	outFileDf.pop(0)
	outName = args.outFile + 'results.tsv'
	outFileDf.fillna(0).to_csv(outName, sep = "\t")

elif mode == "Dir":
	#########################################################################################################
	#### calculate expected probabilities (with direction of change, AE!=EA) and compare with observed ####
	#########################################################################################################

	# Load neutral model
	# Here as prob

	probDf = pd.DataFrame.from_dict(prob, orient="index")
	sum_prob_exp = probDf.sum().sum()

	# From triplets of codons we have to create the file N1_changes_with_direction
	#file_align_list = ["/home/jmontanes/Documents/IQtree_Gene_duplication/Insects/DropboxThings/gain_of_acidic_residues/gain_of_acidic_residues/AlignmentChangesN1_rootSpecies.tsv"]
	fileInDf = pd.read_csv(file_align_list[0], sep="\t", index_col=0).drop("Perf")
	prob_N1 = {}
	for index, row in fileInDf.iterrows():
		if index[0] == index[1] or index[0] == index[2]:
			alDir = ''.join(dict.fromkeys(index))
			try:
				prob_N1[alDir] = prob_N1[alDir] + row["Hits"]
			except:
				prob_N1[alDir] = row["Hits"]

	prob_N1Df = pd.DataFrame.from_dict(prob_N1, orient="index")
	# Remove those that we have 0 probabilities to find it

	prob_N1Df[prob_N1Df.index.isin(probDf[probDf[0] == 0].index)] = 0

	sum_prob_N1 = prob_N1Df.sum().sum()

	scaling_factor = sum_prob_N1 / sum_prob_exp
	outFileDf = probDf.copy()
	fileName = file_align_list[0].split("/")[-1].split(".")[0]
	outFileDf["expected_" + fileName] = outFileDf[0] * scaling_factor
	outFileDf["observed_" + fileName] = prob_N1Df[0]

	firstLetter = outFileDf.index.get_level_values(outFileDf.index.name).str[0].tolist()
	firstType = [aa_type[x] for x in firstLetter]

	secondLetter = outFileDf.index.get_level_values(outFileDf.index.name).str[1].tolist()
	secondType = [aa_type[x] for x in secondLetter]

	Type = []

	for s1, s2 in zip(firstType, secondType):
			Type.append(s1 + "-" + s2)

	outFileDf["Type"] = Type
	outFileDf.pop(0)
	outName = args.outFile + 'results.tsv'
	outFileDf.fillna(0).to_csv(outName, sep = "\t")
else:
	print("Incorrect mode selected, quiting...")
	quit()