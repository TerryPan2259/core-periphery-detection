import csv
import numpy as np
import networkx as nx

from kmalgorithm import km_config, km_modmat

def test1():
	network = nx.karate_club_graph()
	cppairs = km_config(network, significance_level = 0.05)

	assert 'pair_id' in cppairs 
	assert 'core_node' in cppairs 
	assert 'cp_pair_significance' in cppairs 

def test2():
	network = nx.karate_club_graph()

	cppairs = km_modmat(network, significance_level = 0.05)

	assert 'pair_id' in cppairs 
	assert 'core_node' in cppairs 
	assert 'cp_pair_significance' in cppairs 
