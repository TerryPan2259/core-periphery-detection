from abc import ABCMeta, abstractmethod
import abc
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
		""" Private """
		pass

	@abstractmethod
	def _score(self, G, c, x):
		""" Private """
		pass

	def get_pair_id(self):
		"""Get core-periphery pair ID of each node.
	
		
		Returns
		-------
		c : dict
			Key: Node name
			Value: IDs of core-periphery pair to which it belongs. 
			
		"""
		return self.c_
	
	def get_coreness(self):
		"""Get the coreness of each node"""
		return self.x_
		
	def score(*args):
		"""Get score of core-periphery pairs.
		
		Parameters
		----------
		G : Graph object. 
		c : Dict object,  the keys and values of which are the name of node and its ID of belonging core-periphery pair.  
		
		
		Returns
		-------
		q : List. q[i] is the quality of core-periphery pair i.  
			
		"""
		self = args[0]

		if len(args) ==1:
			return self.qs_
		else:
			G = args[1]
			c = args[2]
			x = args[3]
			return self._score(G, c, x)
	
	def _to_edge_list(self, G):
		"""Transform NetworkX object to an edge list. 
		
		Parameters
		----------
		G : Graph object. 
		
		
		Returns
		-------
		node_pairs : (M, 2) numpy array, where M is the number of edges. node_pairs[i,0] and node_pairs[i,1] are the endpoints of the ith edge.  
		w : Mx1 numpy array. w[i] is the weight of the ith edge. 
		node2id : Dict. A function mapping from node name to node id, i.e., node2id[node_name] gives the id.    
		id2node : Dict. A function mapping from node id to node name, i.e., id2node[node_id] gives the node name.
		"""
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
