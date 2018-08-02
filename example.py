import csv
import numpy as np
import pandas as pd
import networkx as nx
import cpalgorithm as cp

G=nx.karate_club_graph()

be = cp.Divisive()
Q = []
be.detect(G)
c = be.get_pair_id()
x = be.is_core()
print(c,x)
