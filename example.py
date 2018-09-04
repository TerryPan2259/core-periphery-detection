import pandas as pd
import networkx as nx
import cpalgorithm as cp

df = pd.read_csv('data/example_edge_list.dat', sep=',') # Load a list of edges (comma-separated file)

G = nx.from_pandas_edgelist(df) # NetworkX graph object

G = nx.karate_club_graph() 
be = cp.MINRES()
be.detect(G)
c = be.get_pair_id()
x = be.get_coreness()

print('Name\tPairID\tCoreness')
for key, value in sorted(c.items(), key=lambda x: x[1]):
	print('%s\t%d\t%f' %(key, c[key], x[key]))
