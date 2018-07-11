import csv
import numpy as np
import pandas as pd
import networkx as nx
import cpalgorithm as cp

G=nx.karate_club_graph()
<<<<<<< HEAD
=======
#G=nx.florentine_families_graph()
#df = pd.read_csv("karate.dat", sep='\t');
#G = nx.from_pandas_edgelist(df, "source", 'target', 'weight')

print(G)
>>>>>>> 5fea595269c406d0e791d280f12d35e3baa3164e
	
be = cp.KM_config()

Q = []
be.detect(G, 1)
c = be.get_pair_id()
x = be.is_core()
<<<<<<< HEAD

print("pair IDs")
print(c)

print("")
print("Role")
print(c)
=======
print(sum(be.score()))

significance, p_values, q_tilde, s_tilde = cp.qstest(c, x, G, be, num_of_thread = 2)
print(c,x)
print(significance, p_values)
>>>>>>> 5fea595269c406d0e791d280f12d35e3baa3164e
