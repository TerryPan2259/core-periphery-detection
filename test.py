import csv
import numpy as np
import networkx as nx

from cpalgorithm import km_config, km_modmat, be, minres

def test1():
	network = nx.karate_club_graph()
	cppairs = km_config(network)

	assert 'pair_id' in cppairs 
	assert 'core_node' in cppairs 

def test2():
	network = nx.karate_club_graph()

	cppairs = km_modmat(network)

	assert 'pair_id' in cppairs 
	assert 'core_node' in cppairs 

def test3():
	network = nx.karate_club_graph()

	cppairs = be(network)

	assert 'pair_id' in cppairs 
	assert 'core_node' in cppairs 

def test4():

	network = nx.karate_club_graph()
	cppairs = minres(network)

	assert 'pair_id' in cppairs 
	assert 'core_node' in cppairs 
