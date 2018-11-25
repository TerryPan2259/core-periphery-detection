import _cpalgorithm as _cp
from .CPAlgorithm import * 

class MINRES(CPAlgorithm):
	"""MINRES algorithm.
	
	MINRES algorithm for finding discrete core-periphery pairs [1], [2]. 
		
	Examples
	--------
	Create this object.

	>>> import cpalgorithm as cpa	
	>>> mrs = cpa.MINRES()
	
	**Core-periphery detection**
	
	Detect core-periphery structure in network G (i.e., NetworkX object):
	
	>>> mrs.detect(G) 
	
	Retrieve the ids of the core-periphery pair to which each node belongs:
	
	>>> pair_id = mrs.get_pair_id() 
	
	Retrieve the coreness:

	>>> coreness = mrs.get_coreness() 

	.. note::

	   This algorithm accepts unweighted and undirected networks only.
	   Also, the algorithm assigns all nodes into the same core-periphery pair by construction, i.e., c[node_name] =0 for all node_name.
	   This algorithm is deterministic, i.e, one obtains the same result at each run. 

	.. rubric:: References

	[1] J. P. Boyd, W. J Fitzgerald, M. C. Mahutga, and D. A. Smith. Computing continuous core/periphery structures for social relations data with MINRES/SVD. Soc.~Netw., 32:125–137, 2010.

	[2] S. Z. W.~ Lip. A fast algorithm for the discrete core/periphery bipartitioning problem. arXiv, pages 1102.5511, 2011.

	"""
	
	def __init__(self):
		self.num_runs = 0 
	
	def detect(self, G):
		"""Detect a single core-periphery pair using the MINRES algorithm.
			
		Parameters
		----------
		G : NetworkX graph object
		
		Examples
		--------

		>>> import networkx as nx
		>>> import cpalgorithm as cpa
		>>> G = nx.karate_club_graph()  # load the karate club network. 
		>>> mrs = cp.MINRES()
		>>> mrs.detect(G)
		
		

		"""

		node_pairs, w, node2id, id2node = self._to_edge_list(G)

		cppairs = _cp.detect_minres(edges=node_pairs, ws=w)
		
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
	
		result = _cp.calc_Q_minres(edges=node_pairs, ws=w, c=_c, x=_x)
		return result[1].tolist()
