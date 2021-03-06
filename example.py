import pandas as pd
import networkx as nx
import cpalgorithm as cp

G = nx.karate_club_graph() 
be = cp.KM_config(num_runs = 100000)
be.detect(G)
c = be.get_pair_id()
x = be.get_coreness()

print('Name\tPairID\tCoreness')
for key, value in sorted(c.items(), key=lambda x: x[1]):
	print('%s\t%d\t%f' %(key, c[key], x[key]))
