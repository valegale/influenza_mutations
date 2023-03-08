import networkx as nx

def compute_single_mutation_support(frequentsets, filtered_nodes):
	
	single_mutation_support = {}
	for index, row in frequentsets.iterrows():
		if len(row["itemsets"]) == 1 and list(row["itemsets"])[0] in filtered_nodes:
			single_mutation_support[list(row["itemsets"])[0]] = round(row["support"], 3)
	return single_mutation_support

def compute_high_support_mutation(frequentsets, filtered_nodes):
	
	high_support = []
	for index, row in frequentsets.iterrows():
		if len(row["itemsets"]) == 1 and list(row["itemsets"])[0] in filtered_nodes and row["support"] > 0.2:
			high_support.append(list(row["itemsets"])[0])
	return high_support

def compute_cooccurring_pair(rules, filtered_nodes, z_metric_threshold = 0.85):
	co_occurring = {}
	for index, row in rules.iterrows():
		if len(row["antecedents"]) == 1 and len(row["consequents"]) == 1:
			if row['zhangs_metric'] > z_metric_threshold and list(row["antecedents"])[0] in filtered_nodes and list(row["consequents"])[0] in filtered_nodes:
				pair = frozenset(row['antecedents']).union(frozenset(row['consequents']))
				if pair not in co_occurring:
					co_occurring[pair] = {}
					co_occurring[pair]["support"] = round(row["support"], 3)
					co_occurring[pair][(list(row["antecedents"])[0] , list(row["consequents"])[0])] = round(row["zhangs_metric"], 3)
				else:
					co_occurring[pair][(list(row["antecedents"])[0] , list(row["consequents"])[0])] = round(row["zhangs_metric"], 3)	
	return co_occurring		

def predict_n_glycosylation_sites(co_occurring, filtered_nodes):	
	data = []
	for pair in co_occurring:
		data.append(tuple(pair))
	
	G = nx.Graph()
	G.add_edges_from(data)
	
	N_sites = set()
	N_sites_reverse = []
	for node in filtered_nodes:
		if node[-1] == "N":		
			for connected_component in nx.connected_components(G):
				if node in connected_component:
					for node_2 in connected_component:
						if node_2[-1] in "ST" and (int(node_2[4:-1]) == (int(node[4:-1]) - 2) or int(node_2[4:-1]) == (int(node[4:-1]) + 2)):
							N_sites.add(frozenset([node]).union(frozenset([node_2])))
							
						else:
							N_sites.add(frozenset([node]))
		if node[3] == "N":		
			N_sites_reverse.append(node)
	#remove subsets
	N_sites_prediction = []
	for mutation in N_sites:
		found = False
		if len(mutation) == 1:
			for pair in N_sites:
				if len(pair) == 2:
					if list(mutation)[0] in pair:
						found = True
			if not found:
				N_sites_prediction.append([list(mutation)[0]])
		else:
			N_sites_prediction.append(list(mutation))
	
	return G, N_sites_prediction, N_sites_reverse
