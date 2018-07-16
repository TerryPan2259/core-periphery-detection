import _cpalgorithm as _cp
from .CPAlgorithm import * 

class SBM(CPAlgorithm):
	"""Stochastic block model.
	
	Core-periphery detection algorithm based on stochastic block models [1]
		
	Examples
	--------
	Create this object.

	>>> import cpalgorithm as cpa	
	>>> sbm = cpa.SBM()
	
	**Core-periphery detection**
	
	Detect core-periphery structure in network G (i.e., NetworkX object):
	
	>>> sbm.detect(G) 
	
	Retrieve the ids of the core-periphery pair to which each node belongs:
	
	>>> pair_id = sbm.get_pair_id() 
	
	Retrieve the coreness:

	>>> coreness = sbm.get_coreness() 
		
	.. note::

	   This algorithm accepts unweighted and undirected networks only.
	   Also, the algorithm assigns all nodes into the same core-periphery pair by construction, i.e., c[node_name] =0 for all node_name.
	   This algorithm is stochastic, i.e., one would obtain different results at each run.

	.. rubric:: Reference

	[1]  X. Zhang, T. Martin, and M. Newman. Identification of core-periphery structure in networks. Phys. Rev. E., 91(3):032803, 2015.

	"""
	
	def __init__(self):
		return
	
	def detect(self, G):
		"""Detect a single core-periphery pair.
	
		Parameters
		----------
		G : NetworkX graph object
		
		Examples
		--------
		>>> import networkx as nx
		>>> import cpalgorithm as cpa
		>>> G = nx.karate_club_graph()  # load the karate club network. 
		>>> sbm = cp.SBM()
		>>> sbm.detect(G)

		"""

		node_pairs, w, node2id, id2node = self._to_edge_list(G)
		
		cppairs = _cp.detect_sbm(edges=node_pairs, ws=w, num_of_runs = 1)
		
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
	
		result = _cp.calc_Q_sbm(edges=node_pairs, ws=w, c=_c, x=_x)

		return result[1].tolist()	
