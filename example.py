import csv
import numpy as np
import pandas as pd
import networkx as nx

#import km_config as kmconfig
import cpalgorithm as cp
#from cpalgorithm import km_config, km_modmat, be

#linkfilename='example_edge_list.txt'
#df = pd.read_csv(linkfilename, sep=' ')
G=nx.karate_club_graph()
#Gi = nx.convert_node_labels_to_integers(G)
node2id = dict(zip(G.nodes, range(len(G.nodes))))
id2node= dict((v,k) for k,v in node2id.items())

nx.relabel_nodes(G, node2id,False)
edges = G.edges(data="weight")	

node_pairs = np.array([ [edge[0], edge[1]] for edge in edges ]).astype(int)
w = np.array([ edge[2] for edge in edges ]).astype(float)

if all(np.isnan(w)):
	nx.set_edge_attributes(G, values =1, name='weight')
	w[:] = 1.0
cppairs = cp.minres(G)

print(cppairs)
