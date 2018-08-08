import _cpalgorithm as _cp
from .CPAlgorithm import * 

class Divisive(CPAlgorithm):
	"""Divisive algorithm.

	An algorithm for finding multiple core-periphery pairs in networks.
	This algorithm partitions a network into communities using the Louvain algorithm. 
	Then, it partitions each community into a core and a periphery using the BE algorithm.
	The quality of a community is computed by that equipped with the BE algorithm.
		
	Parameters
	----------
	num_runs : int
		   Number of runs of the algorithm (optional, default: 10)
		   Run the algorithm num_runs times. Then, this algorithm outputs the result yielding the maximum quality. 
	
	Examples
	--------
	Create this object.

	>>> import cpalgorithm as cpa	
	>>> dv = cpa.Divisive()
	
	**Core-periphery detection**
	
	Detect core-periphery structure in network G (i.e., NetworkX object):
	
	>>> dv.detect(G) 
	
	Retrieve the ids of the core-periphery pair to which each node belongs:
	
	>>> pair_id = dv.get_pair_id() 
	
	Retrieve the coreness:

	>>> coreness = dv.get_coreness() 
		
	.. note::

	   This algorithm accepts unweighted and undirected networks only.
	   This algorithm is stochastic, i.e., one would obtain different results at each run.

	.. rubric:: Reference

        [1] S. Kojaku and N. Masuda. Core-periphery structure requires something else in the network. New Journal of Physics, 20(4):43012, 2018

	"""
	
	def __init__(self, num_runs = 10):
		self.num_runs = num_runs 

	
	def detect(self, G):
		"""Detect a single core-periphery pair using the Borgatti-Everett algorithm.
	
		Parameters
		----------
		G : NetworkX graph object
		
		Examples
		--------
		>>> import networkx as nx
		>>> import cpalgorithm as cpa
		>>> G = nx.karate_club_graph()  # load the karate club network. 
		>>> dv = cpa.Divisive()
		>>> dv.detect(G)

		"""

		node_pairs, w, node2id, id2node = self._to_edge_list(G)
		
		# divide a network into communities
		cppairs = _cp.detect_divisive(edges=node_pairs, ws=w, num_of_runs = self.num_runs)
		
		
		N = len(id2node) 
		self.c_ = dict(zip( [id2node[i] for i in range(N)], cppairs[0].astype(int)))
		self.x_ = dict(zip( [id2node[i] for i in range(N)], cppairs[1]))
		self.Q_ = cppairs[2][0]
		self.qs_ = cppairs[3].tolist()

	def _score(self, G, c, x):
		node_pairs, w, node2id, id2node = self._to_edge_list(G)
	
		N = len(id2node)
		_c = np.array([ c[id2node[i]]  for i in range(N) ])
		_x = np.array([ x[id2node[i]]  for i in range(N) ])
	
		result = _cp.calc_Q_divisive(edges=node_pairs, ws=w, c=_c, x=_x)

		return result[1].tolist()	
