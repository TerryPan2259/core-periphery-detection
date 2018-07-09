import csv
import numpy as np
import pandas as pd
import networkx as nx
import cpalgorithm as cp

G=nx.karate_club_graph()
#G=nx.florentine_families_graph()
	
be = cp.KM_config()

Q = []
for i in range(1000):
	be.detect(G)
	c = be.get_pair_id()
	x = be.is_core()
	Q.append(sum(be.score()))

print(sum(Q) / 100)
#significance, p_values, q_tilde, s_tilde = cp.qstest(c, x, G, be)

#print(significance, p_values)
