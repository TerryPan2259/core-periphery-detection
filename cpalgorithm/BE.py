import _cpalgorithm as _cp
from .CPAlgorithm import * 

class BE(CPAlgorithm):
	"""Borgatti Everett algorithm.

	An algorithm for finding single core-periphery pair in networks.
		
	Parameters
	----------
	num_runs : int
		   Number of runs of the algorithm (optional, default: 10)
		   Run the algorithm num_runs times. Then, this algorithm outputs the result yielding the maximum quality. 
	
	Examples
	--------
	Create this object.

	>>> import cpalgorithm as cpa	
	>>> be = cpa.BE()
	
	**Core-periphery detection**
	
	Detect core-periphery structure in network G (i.e., NetworkX object):
	
	>>> be.detect(G) 
	
	Retrieve the ids of the core-periphery pair to which each node belongs:
	
	>>> pair_id = be.get_pair_id() 
	
	Retrieve the coreness:

	>>> coreness = be.get_coreness() 
		
	.. note::

	   This algorithm accepts unweighted and undirected networks only.
	   Also, the algorithm assigns all nodes into the same core-periphery pair by construction, i.e., c[node_name] =0 for all node_name.
	   This algorithm is stochastic, i.e., one would obtain different results at each run.

	.. rubric:: Reference

	[1] S. P. Borgatti and M. G. Everett. Models of core/periphery structures. Soc.~Netw., 21(4):375–395, 2000.

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
		>>> be = cp.BE()
		>>> be.detect(G)

		"""

		node_pairs, w, node2id, id2node = self._to_edge_list(G)

		cppairs = _cp.detect_be(edges=node_pairs, ws=w, num_of_runs = self.num_runs)
		
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
	
		result = _cp.calc_Q_be(edges=node_pairs, ws=w, c=_c, x=_x)

		return result[1].tolist()	
