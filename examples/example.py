"""
==========================
Borgatti Everett algorithm
==========================
Example of parallel implementation of betweenness centrality using the
multiprocessing module from Python Standard Library.
The function betweenness centrality accepts a bunch of nodes and computes
the contribution of those nodes to the betweenness centrality of the whole
network. Here we divide the network in chunks of nodes and we compute their
contribution to the betweenness centrality of the whole network.
This doesn't work in python2.7.13. It does work in 3.6, 3.5, 3.4, and 3.3.
It may be related to this:
https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map
"""


import csv
import numpy as np
import pandas as pd
import networkx as nx
import _cpalgorithm as _cp
import cpalgorithm as cp

G=nx.karate_club_graph()
#G=nx.florentine_families_graph()
#df = pd.read_csv("karate.dat", sep='\t');
#G = nx.from_pandas_edgelist(df, "source", 'target', 'weight')

be = cp.BE()

Q = []
be.detect(G)
c = be.get_pair_id()
x = be.is_core()

print(sum(be.score()))

#significance, p_values, q_tilde, s_tilde = cp.qstest(c, x, G, be, num_of_thread = 4, null_model = cp.erdos_renyi)
print(c,x)
#print(significance, p_values)
