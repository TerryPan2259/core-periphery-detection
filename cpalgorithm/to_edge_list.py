
def to_edge_list(G):
	node2id = dict(zip(G.nodes, range(len(G.nodes))))
	id2node= dict((v,k) for k,v in node2id.items())

	nx.relabel_nodes(G, node2id,False)
	edges = G.edges(data="weight")	

	node_pairs = np.array([ [edge[0], edge[1]] for edge in edges ]).astype(int)
	w = np.array([ edge[2] for edge in edges ]).astype(float)

	if all(np.isnan(w)):
		nx.set_edge_attributes(G, values =1, name='weight')
		w[:] = 1.0
	
	nx.relabel_nodes(G,id2node,False)

	return node_pairs, w, node2id, id2node
