from pyvis.network import Network
import networkx as nx
from antigenic_sites import antigenic_sites


def compute_filtered_nodes(frequentsets, recent_mutations):
	
	filtered_nodes = {}
	reverse_mutations = []
	for index, row in frequentsets.iterrows():
		if len(row["itemsets"]) == 1:
			mutation = list(row["itemsets"])[0]
			X = mutation[3]
			Y = mutation[-1]
			reverse_mutation = mutation[0:3] + Y + mutation[4:len(mutation)-1] + X
			if reverse_mutation in recent_mutations:
				if mutation[3] not in "XB" and mutation[-1] not in "XB":
					reverse_mutations.append(reverse_mutation)
			elif reverse_mutation not in filtered_nodes:
				filtered_nodes[mutation] = row["support"]
			else: 
				if row["support"] > filtered_nodes[reverse_mutation]:
					filtered_nodes[mutation] = row["support"]
					del(filtered_nodes[reverse_mutation])
	return filtered_nodes, reverse_mutations
	
def network_complete(frequentsets, rules, recent_mutations, metric = 'support', z_metric_threshold = 0.85):
	
	filtered_nodes, disappearing_mutations = compute_filtered_nodes(frequentsets, recent_mutations)
					
	# Creating a network for the analysis (networkx) and a network for visualization (pyvis)
	net = Network("1000px", "2000px", directed=True, bgcolor = "white", select_menu = False, filter_menu=False)
	net.repulsion()
	G = nx.DiGraph()
	
	# adding nodes 
	for protein_mutation in filtered_nodes:
		protein = protein_mutation[0:2]
		mutation = protein_mutation[:3] + protein_mutation[4:-1]
		
		if protein == "ha":
			group = 1
		else:
			group = 2
		
		if mutation in antigenic_sites:
			border_size = 4
		else:
			border_size = 0
					
		net.add_node(protein_mutation, title = "Support: {}".format(round(filtered_nodes[protein_mutation],3)) , group = group, borderWidth = border_size)
		
	#adding edges from rules
	
	i = 1
	for index, row in rules.iterrows():
			
		X = row['antecedents']
		Y = row['consequents']
		all_real = True
		for protein_mutation in X:
			if protein_mutation not in filtered_nodes:
				all_real = False
		for protein_mutation in Y:
			if protein_mutation not in filtered_nodes:
				all_real = False
		if all_real:
			net.add_node(str(i), shape = 'box', title = "Support: {}\nConfidence: {}\nZhang's metric: {}".format(round(row['support'],3), round(row['confidence'],3), round(row['zhangs_metric'],3)), color = "f0d1ef")
			G.add_node(str(i))
			
			if row['zhangs_metric'] > z_metric_threshold:
				color = "6CC776"	
			
			else:
				color = "grey"
			for protein_mutation in X:
				G.add_node(protein_mutation)
				G.add_edge(protein_mutation, str(i))
				net.add_edge(protein_mutation, str(i), value = row[metric], color = color)
			
			for protein_mutation in Y:
				G.add_node(protein_mutation)
				G.add_edge(str(i), protein_mutation)
				net.add_edge(str(i), protein_mutation, value = row[metric], color = color)
			
			i += 1
	
	return net, G, filtered_nodes, disappearing_mutations

def in_degree_analysis(G):
	# analysis out degree vs in degree, nodes with large out_degree are more likely to be 
	# key mutations
	difference_in_out = {}
	for (mutation, degree) in G.in_degree:
		if mutation[0] in ('hn'): #excluding rules-nodes
			difference_in_out[mutation] = degree
	
	for (mutation, degree) in G.out_degree:
		if mutation[0] in ('hn'): #excluding rules-nodes
			difference_in_out[mutation] -= degree
			
	x = {key:value for (key, value) in difference_in_out.items() if value >= 2}
	
	return (sorted(x.keys(), key=lambda item: item[0], reverse= True), x)

