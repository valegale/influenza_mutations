from transactions import generate_transaction_list
from frequent_itemsets import frequent_itemsets_fpgrowth
from association_rules import generate_association_rules
from network_analysis import network_complete, in_degree_analysis
from extract_cooccurring import compute_single_mutation_support, compute_cooccurring_pair, predict_n_glycosylation_sites
from create_yattag import create_document, create_clusters_summary_document
import networkx as nx
import random
import os


#### Parameters
bin_start = "18_19"
bin_end = "19_20"
path = "./../H3N2/removed_duplicates"

N = 250
max_size_itemset = 2
min_support = 0.05
min_confidence = 0.5
random_seed = 100

only_antigenic_sites = False
metric = "support"
threshold = 15
####

## Generate transactions
if random_seed:
	random.seed(random_seed)
	all_transactions = generate_transaction_list(path, bin_start, bin_end, N, threshold)
else:
	all_transactions = generate_transaction_list(path, bin_start, bin_end, N, threshold)

def padding_zero(i):
	return str(i).zfill(2)

# create a result folder if it doesn't exist
path_results = "results_temp"
if not os.path.exists(path_results):
   os.makedirs(path_results)
   print("A new folder called {} has been created".format(path_results))
   
# looping through each pair of flu season
cluster_dictionary = []
recent_mutations = {}
j = 0
for i in range(int(bin_start[0:2]), int(bin_end[0:2])):
	bin1 = "{}_{}".format(padding_zero(i), padding_zero(i+1))
	bin2 = "{}_{}".format(padding_zero(i+1), padding_zero(i+2))
	transactions = all_transactions[(j*N):((j+1)*N)] 
	j += 1
	
	print ("From flu season {} to flu season {}".format(bin1, bin2))
	## Run fpgrowth to compute the frequent itemsets
	frequentsets = frequent_itemsets_fpgrowth(transactions, min_supp=min_support, max_len = max_size_itemset)
	
	print("Number of frequent itemsets: {}".format(str(len(frequentsets))))
	
	## Generate association rules
	rules = generate_association_rules(frequentsets, metric="confidence", min_threshold=min_confidence, only_antigenic_sites = only_antigenic_sites, only_transition_sites = False)

	print("Number of rules: {}".format(str(len(rules))))
	
	net, G, filtered_nodes, reverse_mutations = network_complete(frequentsets, rules, recent_mutations)

	recent_mutations = filtered_nodes  # keeping track of the mutation of the previous year for the next flu season
	key_mutations, key_mutations_dictionary = in_degree_analysis(G) #key mutations: indegree - outdgree > 2
	single_mutation_support = compute_single_mutation_support(frequentsets, filtered_nodes)
	co_occurring = compute_cooccurring_pair(rules, filtered_nodes, z_metric_threshold = 0.85)
	G_cooccurring, N_sites_prediction, N_sites_reverse = predict_n_glycosylation_sites(co_occurring, filtered_nodes)

	
	for connected_component in nx.connected_components(G_cooccurring):
		support_dictionary = {}
		for mutation in connected_component:
			support_dictionary[mutation] = single_mutation_support[mutation]
		cluster_dictionary.append(([bin1 + "-" + bin2], support_dictionary))
	

	name_file = "/{}_{}.html".format(bin1, bin2)
	create_document(name_file, path_results, bin1, bin2, min_support, net, G_cooccurring, key_mutations, single_mutation_support, co_occurring, N_sites_prediction, N_sites_reverse)
	

#formatting for cluster table
data_table_clusters = []
for flu_season, clusters in cluster_dictionary:		
	
	mutations = []
	mutations_support = []
	only_position_ha = []
	only_position_na = []
	for mutation in clusters:
		mutations.append(mutation)
		mutations_support.append((clusters[mutation]))
		if mutation[:2] == "ha":
			only_position_ha.append(mutation[4:-1])
		if mutation[:2] == "na":
			only_position_na.append(mutation[4:-1])
	data_table_clusters.append((flu_season, mutations, mutations_support, only_position_ha, only_position_na))


create_clusters_summary_document(path_results, data_table_clusters)