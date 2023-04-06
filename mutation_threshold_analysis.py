from Bio import SeqIO

path = "./../H1N1/removed_duplicates"
bin1 = "14_15"
bin2 = "15_16"

def find_mutations(seq1, seq2, protein = "ha"):
	return [protein + "_" + seq1[i - 1] + str(i) + seq2[i - 1] for i in range(1, len(seq1) + 1) if seq1[i - 1] != seq2[i - 1] and seq1[i - 1] != "-" and seq2[i - 1] != "-"]

ha_protein_file_bin1 = "{}/ha_proteins/ha_{}.fa".format(path, bin1)
na_protein_file_bin1 = "{}/na_proteins/na_{}.fa".format(path, bin1)
	

sequence1 = "A/Minnesota/32/2015"
sequence2 = "A/Nicaragua/6670_10/2015"

for seq_record in SeqIO.parse(ha_protein_file_bin1, "fasta"):
	if seq_record.id == sequence1:
		seq1_ha = seq_record.seq

for seq_record in SeqIO.parse(na_protein_file_bin1, "fasta"):
	if seq_record.id == sequence1:
		seq1_na = seq_record.seq

ha_protein_file_bin2 = "{}/ha_proteins/ha_{}.fa".format(path, bin2)
na_protein_file_bin2 = "{}/na_proteins/na_{}.fa".format(path, bin2)

for seq_record in SeqIO.parse(ha_protein_file_bin2, "fasta"):
	if seq_record.id == sequence2:
		seq2_ha = seq_record.seq

for seq_record in SeqIO.parse(na_protein_file_bin2, "fasta"):
	if seq_record.id == sequence2:
		seq2_na = seq_record.seq
		

mutations_ha = find_mutations(seq1_ha, seq2_ha)
print (mutations_ha)

mutations_na = find_mutations(seq1_na, seq2_na)
print (mutations_na)

print (len(mutations_ha) + len(mutations_na))
