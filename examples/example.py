"""Detect core-periphery structure in empirical networks.

"""
import cpalgorithm as cp
import scipy.stats as stats
import networkx as nx
import cpalgorithm as cpa
import numpy as np 
import sys

def profiling_cp_structure(G, alg, null_model):
	
	# cp detection
	alg.detect(G) 

	# retrieve result	
	pair_id = alg.get_pair_id()
	coreness = alg.get_coreness()

	# statistical test	
	pair_id, coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, alg, null_model=null_model, num_of_rand_net = 500 )
	print(p_values)
	# Number of pairs
	C = len(significance)
	
	# number of nodes
	N = nx.number_of_nodes(G)

	# filter out all None values	
	c = [v for k,v in pair_id.items() if v != None]
	x = [v for k,v in coreness.items() if v != None]
	if len(c) == 0:
		Csig = 0
		fc = 0
		x = 0
		corr = 0	
	else:
		# Number of significant pairs
		Csig = max(c)
		
		# Average coreness value (in case of disc. cp, this is the fraction of core nodes)
		fc = sum(x) / float(N)
	
		# compute the correlation between degree and the coreness 
		deg = [G.degree(k) for k,v in coreness.items() if v != None]
		cness = [v for k,v in coreness.items() if v != None]
		corr = stats.kendalltau(deg, cness).correlation
	
	return C, Csig, fc, corr


# loading graphs
Glist = {\
	"karate":nx.karate_club_graph(),\
	"dolphins":nx.read_gml("data/dolphins.gml"),\
	"book":nx.read_gml("data/polbooks.gml"),\
	"celegans":nx.read_gml("data/celegansneural.gml"),\
	#"jazz":nx.read_pajek("data/jazz.net"),\
	"lesmiserables":nx.read_gml("data/lesmiserables.gml"),\
	#"polblogs":nx.read_gml("data/polblogs.gml"),\
	#"power":nx.read_gml("data/power.gml", label='id')\
	};

alg = cpa.BE()
for name, G in Glist.items():
	print(name)
	G = G.to_undirected()
	print(profiling_cp_structure(G, alg, cpa.erdos_renyi))
