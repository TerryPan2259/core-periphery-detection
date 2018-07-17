import csv
import numpy as np
import pandas as pd
import networkx as nx
import cpalgorithm as cp
from scipy.sparse import diags
from scipy.sparse.linalg import eigs

G=nx.karate_club_graph()
N = nx.number_of_nodes(G)
M = nx.number_of_edges(G) 
A = nx.to_scipy_sparse_matrix(G)
deg = np.array([d[1] for d in G.degree()])
denom =  1.0 / deg
denom[np.isnan(denom)] = 0
T = diags(denom) * A

d, v = eigs(T, k=1, sigma=1e-2);

iscore = np.sign(v)>0
print(iscore)
