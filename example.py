import csv
import numpy as np
import pandas as pd
import networkx as nx
import cpalgorithm as cp

G=nx.karate_club_graph()

be = cp.Rombach()
Q = []
be.detect(G)
c = be.get_pair_id()
x = be.is_core()
print(sum(be.score()))
print(c,x)

be = cp.LapSgnCore()
Q = []
be.detect(G)
c = be.get_pair_id()
x = be.is_core()
print(sum(be.score()))
print(c,x)
#significance, p_values, q_tilde, s_tilde = cp.qstest(c, x, G, be, num_of_thread = 1, null_model = cp.erdos_renyi)
#print(significance, p_values)
