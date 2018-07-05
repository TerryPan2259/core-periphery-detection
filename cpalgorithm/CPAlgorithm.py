from abc import ABCMeta, abstractmethod
import networkx as nx
import numpy as np

class CPAlgorithm(metaclass=ABCMeta):

	def __init__(self):
		self.x_ = [] 
		self.c_ = []
		self.Q_ = []
		self.qs_ = [] 

	@abstractmethod
	def detect(self):
		pass


	@abstractmethod
	def _score(self, G, c, x):
		pass

	@abstractmethod
	def significance(self):
		pass

	def get_pair_id(self):
		return self.c_
	
	def is_core(self):
		return self.x_
		
	def score(*args):
		self = args[0]

		if len(args) ==1:
			return self.Q_, self.qs_
		else:
			G = args[1]
			c = args[2]
			x = args[3]
			return self._score(G, c, x)
	
	def to_edge_list(self, G):
		node2id = dict(zip(G.nodes, range(len(G.nodes))))
		id2node= dict((v,k) for k,v in node2id.items())
	
		nx.relabel_nodes(G, node2id,False)
		edges = G.edges(data="weight")	
	
		node_pairs = np.array([ [edge[0], edge[1]] for edge in edges ]).astype(int)
		w = np.array([ edge[2] for edge in edges ]).astype(float)
	
		if all(np.isnan(w)):
			nx.set_edge_attributes(G, values =1, name='weight')
			w[:] = 1.0
		
		nx.relabel_nodes(G,id2node,False)
	
		return node_pairs, w, node2id, id2node
