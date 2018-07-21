import _cpalgorithm as _cp
from .CPAlgorithm import * 

class KM_config(CPAlgorithm):
	"""Kojaku-Masuda algorithm with the configuration model.
	
	This algorithm finds multiple core-periphery pairs in networks. 
	In the detection of core-periphery pairs, the configuration model is used as the null model. 	
	
	Parameters
	----------
	num_runs : int
		Number of runs of the algorithm  (optional, default: 1).  

	Examples
	--------
	Create this object.

	>>> import cpalgorithm as cpa	
	>>> km = cpa.KM_config()
	
	**Core-periphery detection**
	
	Detect core-periphery structure in network G (i.e., NetworkX object):
	
	>>> km.detect(G) 
	
	Retrieve the ids of the core-periphery pair to which each node belongs:
	
	>>> pair_id = km.get_pair_id() 
	
	Retrieve the coreness:

	>>> coreness = km.get_coreness() 
		
	.. note::

	   This algorithm can accept unweighted and weighted networks.

	.. rubric:: Reference

        [1] S. Kojaku and N. Masuda. Core-periphery structure requires something else in the network. New Journal of Physics, 20(4):43012, 2018

	"""
	
	
	def __init__(self, num_runs = 10):
		self.num_runs = num_runs 
	
	def detect(self, G):
		"""Detect multiple core-periphery pairs.
	
		Parameters
		----------
		G : NetworkX graph object
		
		Examples
		--------
		>>> import networkx as nx
		>>> import cpalgorithm as cpa
		>>> G = nx.karate_club_graph()  # load the karate club network. 
		>>> km = cp.KM_config() # label switching algorithm
		>>> km.detect(G)

		"""

		node_pairs, w, node2id, id2node = self._to_edge_list(G)

		cppairs = _cp.detect_config(edges=node_pairs, ws=w, num_of_runs = self.num_runs)

		N = len(node2id)	
		self.c_ = dict(zip( [id2node[i] for i in range(N)], cppairs[0].astype(int)))
		self.x_ = dict(zip( [id2node[i] for i in range(N)], cppairs[1]))
		self.Q_ = cppairs[2][0]
		self.qs_ = cppairs[3].tolist()

	
	def _score(self, G, c, x):

		node_pairs, w, node2id, id2node = self._to_edge_list(G)
	
		N = len(id2node)
		_c = np.array([ c[id2node[i]]  for i in range(N) ])
		_x = np.array([ x[id2node[i]]  for i in range(N) ])
	
		result = _cp.calc_Q_config(edges=node_pairs, ws=w, c =_c, x = _x)

		return result[1].tolist()
	
	def significance(self):
		return self.pvalues	
