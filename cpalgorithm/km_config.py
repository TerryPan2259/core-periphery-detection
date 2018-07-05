import _cpalgorithm as cp 
import numpy as np
import networkx as nx
import qstest

def km_config(G, num_of_runs = 10, significance_level = 0.05):
	if isinstance(G,nx.classes.graph.Graph) == False:
		print("Pass Networkx.classes.graph.Graph object")
		return

	node_pairs, w, node2id, id2node = to_edge_list(G)

	cppairs = cp.detect_config(edges=node_pairs, ws=w, num_of_runs = num_of_runs)

	if significance_level < 1.0
	
		for cid in range(max(cppiars[0])):
			sl = cid==cppairs[0]
				
	
		def my_qfunc(network, nodes):
			cp.eval_config(edges = node_pairs, ws=w, cppairs[0])
    			return network.subgraph(nodes).size()
		
		qstest.qstest(G, cppairs,  	
		
		
	

	return {'pair_id': dict(zip( range(len(cppairs[0])), cppairs[0].astype(int))),\
		'core_node': dict(zip( range(len(cppairs[0])), cppairs[1].astype(bool)))}

	#return {'pair_id':, 'is_core_node':cppairs[1], 'cp_pair_significance':cppairs[2]}	
