from Bio import SeqIO
import random

def generate_transactions(path, bin1, bin2, N, threshold = False, aminoacid_name = True):
	# parameters
	# path: the path to a folder that includes two folders with all the 
	#       amino acid sequences for ha and na divided into flu season bins, ha_proteins and na_proteins. 
	# bin1, bin2: name of the bins that are compared (example 08_09 and 09_10)
	# N:    number of times the sequences are extracted randomly with replacement from the bins

	ha_protein_file_bin1 = "{}/ha_proteins/ha_{}.fa".format(path, bin1)
	na_protein_file_bin1 = "{}/na_proteins/na_{}.fa".format(path, bin1)
	
	sequences_bin_1 = {}
	for seq_record in SeqIO.parse(ha_protein_file_bin1, "fasta"):
		sequences_bin_1[seq_record.id] = {}
		sequences_bin_1[seq_record.id]["ha_sequence"] = seq_record.seq
	
	for seq_record in SeqIO.parse(na_protein_file_bin1, "fasta"):
		sequences_bin_1[seq_record.id]["na_sequence"] = seq_record.seq
	
	ha_protein_file_bin2 = "{}/ha_proteins/ha_{}.fa".format(path, bin2)
	na_protein_file_bin2 = "{}/na_proteins/na_{}.fa".format(path, bin2)
	
	sequences_bin_2 = {}
	for seq_record in SeqIO.parse(ha_protein_file_bin2, "fasta"):
		sequences_bin_2[seq_record.id] = {}
		sequences_bin_2[seq_record.id]["ha_sequence"] = seq_record.seq
	
	for seq_record in SeqIO.parse(na_protein_file_bin2, "fasta"):
		sequences_bin_2[seq_record.id]["na_sequence"] = seq_record.seq

	
	def find_mutations(seq1, seq2, aminoacid_name, protein = "ha"):
		if aminoacid_name:
			return [protein + "_" + seq1[i - 1] + str(i) + seq2[i - 1] for i in range(1, len(seq1) + 1) if seq1[i - 1] != seq2[i - 1] and seq1[i - 1] != "-" and seq2[i - 1] != "-"]
		else:
			return [protein + "_" + str(i) for i in range(1, len(seq1) + 1) if seq1[i - 1] != seq2[i - 1] and seq1[i - 1] != "-" and seq2[i - 1] != "-"]

			
	list_transactions = []
	for i in range(N):
		if threshold:
			number_mutations = threshold
			
			while threshold <= number_mutations:
				
				sequence1 = random.choice(list(sequences_bin_1.keys()))
				sequence2 = random.choice(list(sequences_bin_2.keys()))
				mutations_ha = find_mutations(sequences_bin_1[sequence1]["ha_sequence"], sequences_bin_2[sequence2]["ha_sequence"], aminoacid_name)
				mutations_na = find_mutations(sequences_bin_1[sequence1]["na_sequence"], sequences_bin_2[sequence2]["na_sequence"], aminoacid_name, "na")
				number_mutations = len(mutations_ha + mutations_na)
			
			list_transactions.append(mutations_ha + mutations_na)
			
			
		else:
			sequence1 = random.choice(list(sequences_bin_1.keys()))
			sequence2 = random.choice(list(sequences_bin_2.keys()))
			mutations_ha = find_mutations(sequences_bin_1[sequence1]["ha_sequence"], sequences_bin_2[sequence2]["ha_sequence"], aminoacid_name)
			mutations_na = find_mutations(sequences_bin_1[sequence1]["na_sequence"], sequences_bin_2[sequence2]["na_sequence"], aminoacid_name, "na")

			list_transactions.append(mutations_ha + mutations_na)
				
	return list_transactions

def generate_transaction_list(path, first_bin, last_bin, N, threshold):
	# generate transactions from all consecutive bins starting from first_bin until last_bin
	# example first_bin = 08_09, last_bin = 13_14
	
	def padding_zero(i):
		return str(i).zfill(2)
		
	transactions = []
	for i in range(int(first_bin[0:2]), int(last_bin[0:2])):
		bin1 = "{}_{}".format(padding_zero(i), padding_zero(i+1))
		bin2 = "{}_{}".format(padding_zero(i+1), padding_zero(i+2))
		new_list_transactions = generate_transactions(path, bin1, bin2, N, threshold)
		transactions = transactions + new_list_transactions
	
	return transactions
	
def generate_transaction_list_only_position(path, first_bin, last_bin, N, threshold):
	# generate transactions from all consecutive bins starting from first_bin until last_bin
	# example first_bin = 08_09, last_bin = 13_14. This function only return the position and 
	# NOT which aminoacid is mutating. Currently not used in the tool.
	
	aminoacid_name = False
	def padding_zero(i):
		return str(i).zfill(2)
		
	transactions = []
	for i in range(int(first_bin[0:2]), int(last_bin[0:2])):
		bin1 = "{}_{}".format(padding_zero(i), padding_zero(i+1))
		bin2 = "{}_{}".format(padding_zero(i+1), padding_zero(i+2))
		new_list_transactions = generate_transactions(path, bin1, bin2, N, threshold, aminoacid_name)
		transactions = transactions + new_list_transactions
	
	return transactions
	