import networkx as nx 

G = nx.read_gml('polblogs.gml')
H = nx.convert_node_labels_to_integers(G)
nx.write_edgelist(H,'polblog.txt')
print(nx.info(H))
