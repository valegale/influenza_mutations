from antigenic_sites import antigenic_sites
import networkx as nx
from yattag import Doc

def create_document(name_file, bin_1, bin_2, min_support, net, G, key_mutation, single_mutation_support, co_occurring, N_sites_prediction, N_sites_reverse):
	doc, tag, text, line = Doc().ttl()
	with tag('head'):
		with tag('style'):
			doc.asis('p {color:#B5301E; display:inline}')
			
	with tag('html'):
	    
		
		with tag('h1'):
			text('Co-occurring mutations in influenza virus')
		text('Analysis of the mutations occurring from flu season {} to flu season {}'.format(bin_1.replace("_", "/" ), bin_2.replace("_", "/" )))
		with tag('h2'):
			text('Co-occurring mutation network')
			
		doc.asis(net.generate_html())
		text('The network includes all mutations with ')
		with tag('b'):
			text('minimum support >=' + str(min_support))
		text(" and ")
		with tag('b'):
			text('minimum confidence >= 0.5. ')		
		text("Each pink box represents a rule. Support, Confidence, and Zhang's metric can be visualized by pointing the cursor directly on the node. ")
		doc.stag('br')
		text('Each mutation is a node and the support for the individual mutation can be visualized by pointing the cursor on the node. HA and NA mutation nodes have different colors. ')
		doc.stag('br')
		text('Mutations in the antigenic regions have a colored border around the node. ')
		doc.stag('br')
		text("Each rule has directed edges coming out and in, according to the direction of the rule. The thickness of the edge indicates the support of the rule and rules with a high Zhang's metric (likely to be associated) are indicated with a green line.")
			
		
		with tag('h1'):
			text('Analysis and prediction')
		
		text('Analysis of the most relevant mutations. Mutations in red occurred in the antigenic region.')
		
		with tag('h2'):
			text('Possible drivers of evolution')
		text('This list is obtained by extracting the nodes with indegree - outdegree >= 2. ')
		text('Mutations with a higher out degree are likely to be key mutations of the virus, as this means that they cannot be predicted by the occurrence of other mutations, but many other mutations are associated with them.')
		# insert predicted drivers
		with tag('table', klass='table table-bordered table-responsive table-striped'):
			with tag('thead',  klass='thead-light'):
				with tag('tr'):
					with tag('th'):
						text('Protein') 
					with tag('th'):
						text('Mutation') 
					with tag('th'):
						text('Aminoacid position') 
					with tag('th'):
						text('Support') 
						
				for mutation in key_mutation:
					
					with tag('tr'):
						with tag('td'):
							if mutation[:2] == "ha":
								text('Hemagglutinin')
							if mutation[:2] == "na":
								text('Neuraminidase')
						with tag('td'):
							text(mutation[3] + " -> "  + mutation[-1]) 
						with tag('td'):
							if mutation[:3] + mutation[4:-1] in antigenic_sites:
								line('p', mutation[4:-1] )
							else:
								text(mutation[4:-1]) 
						with tag('td'):
							text(single_mutation_support[mutation]) 
		with tag('h2'):
			text('List of predicted co-occurring mutations')
		text("The following list includes all pair of association that has a high Zhang's metric")
		with tag('table', klass='table table-bordered table-responsive table-striped'):
			with tag('thead',  klass='thead-light'):
				with tag('tr'):
					with tag('th'):
						text('Mutations')
					with tag('th'):
						text('Support') 
					with tag('th'):
						text('Zhangs metrics') 
	
				
				while len(co_occurring) > 0:
					highest = max(co_occurring, key=lambda key: co_occurring[key]['support'])
					values = co_occurring.pop(highest)				
					support = values.pop('support')
					
					z_metric = []
					for pair in values:
						z_metric.append(values[pair])
					with tag('tr'):
						with tag('td'):
							if len(z_metric) == 2:
								if pair[1][:3] + pair[1][4:-1] in antigenic_sites:
									
									line('p', pair[1])
									text("  <-> ")
								else: 
									text(pair[1] + "  <-> ")
								if pair[0][:3] + pair[0][4:-1] in antigenic_sites:
									line('p', pair[0])
								else:
									text(pair[0])
								#text(pair[1] + "  <-> " + pair[0])
							if len(z_metric) == 1:
								text(pair[0] + "  -> " + pair[1])
						with tag('td'):
							text(support) 
						with tag('td'):
							text(str(z_metric) )	
		
	
		with tag('h3'):
			text('Co-occurring clusters of mutations')
		text('Connected component network analysis is used to identify the subgraphs in which each pair of nodes is connected with each other via a path.')
			
		with tag('table', klass='table table-bordered table-responsive table-striped'):
			with tag('thead',  klass='thead-light'):
				with tag('tr'):
					with tag('th'):
						text('Clusters of mutations')
	
	
				for connected_component in nx.connected_components(G):
					with tag('tr'):
						with tag('td'):
							
							for mutation in connected_component:
								
								if mutation[:3] + mutation[4:-1] in antigenic_sites:
									line('p', mutation)
								else:
									text(mutation)
								text(', ')
								
											
		with tag('h2'):
			text('N-glycosylation site prediction')
			with tag('table', klass='table table-bordered table-responsive table-striped'):
				with tag('thead',  klass='thead-light'):
					with tag('tr'):
						with tag('th'):
							text('Emerging N-glycosylation sites')
						with tag('th'):
							text('Support')
					
				
					for n_site in N_sites_prediction:
						with tag('tr'):
							with tag('td'):
								text(str(n_site))
							with tag('td'):
								list_support = [single_mutation_support[mutation] for mutation in n_site]
								text(str(list_support))	
	
	
			with tag('table', klass='table table-bordered table-responsive table-striped'):
				with tag('thead',  klass='thead-light'):
					with tag('tr'):
						with tag('th'):
							text('Disappearing N-glycosylation sites')
						with tag('th'):
							text('Support')
	
	
					for n_site in N_sites_reverse:
						with tag('tr'):
							with tag('td'):
								text(str(n_site))					
							with tag('td'):
								text(single_mutation_support[n_site])	
								
									
									
	
	f = open("results_html/" + name_file, "w")
	f.write(doc.getvalue())
	f.close()

	
def create_clusters_summary_document(data_table_clusters, name_file = "co_occurring_clusters.html"):
	doc, tag, text = Doc().tagtext()
	
	#
	# create the CSS style
	style = '''
	    table {
	        border-collapse: collapse;
	        font-size: 1em;
	    }
	    th, td {
	        border: 1px solid #ccc;
	        text-align: center;
	    }
	    th {
	        background-color: #eee;
	        font-weight: bold;
	    }
	    header {
	        text-align: left;
	        font-size: 1.5em;
	        font-weight: bold;
	        margin-bottom: 1em;
	    }
	'''
	
	# add the style to the document
	doc.asis('<!DOCTYPE html>')
	with tag('html'):
		with tag('head'):
			with tag('title'):
				text('Clusters results')
			with tag('style'):
				text(style)
			with tag('body'):
				with tag('table'):
					with tag('tr'):
						with tag('th'):
							text('Flu seasons')
						with tag('th'):
							text('Clusters of mutations')
						with tag('th'):
							text('Support')
						with tag('th'):
							text('HA position')
						with tag('th'):
							text('NA position')
					doc.stag('thead')
					with tag('tbody'):
		                # add the table data rows
						for row in data_table_clusters:
							with tag('tr'):
								with tag('td'):
									for obj in row[0]:
										text(str(obj) + " ")
								with tag('td'):
									for obj in row[1]:
										text(str(obj) + " ")
								with tag('td'):
									for obj in row[2]:
										text(str(obj) + " ")
								with tag('td'):
									for obj in row[3]:
										text(str(obj) + " ")
								with tag('td'):
									for obj in row[4]:
										text(str(obj) + " ")
		        # create the download button
				with tag('button', onclick='downloadTable()'):
					text('Download as CSV')
		        # add the JavaScript function to download the table as CSV
				doc.asis('''
		            <script>
		            function downloadTable() {
		                // get the table data as a CSV string
		                var csv = 'Flu_season,mutations,mutation_support, ha_position, na_position\\n';
		                var rows = document.getElementsByTagName('tr');
		                for (var i = 1; i < rows.length; i++) {
		                    var cells = rows[i].getElementsByTagName('td');
		                    for (var j = 0; j < cells.length; j++) {
		                        csv += cells[j].textContent;
		                        if (j < cells.length - 1) {
		                            csv += ',';
		                        }
		                    }
		                    csv += '\\n';
		                }
		                // create a link to download the CSV file
		                var link = document.createElement('a');
		                link.setAttribute('href', 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv));
		                link.setAttribute('download', 'table.csv');
		                link.style.display = 'none';
		                document.body.appendChild(link);
		                link.click();
		                document.body.removeChild(link);
		            }
		            </script>
		        ''')
		
			

									
	f = open("results_html/" + name_file, "w")
	f.write(doc.getvalue())
	f.close()
	

	f.close()
	