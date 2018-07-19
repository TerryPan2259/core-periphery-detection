import _cpalgorithm as _cp
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from .CPAlgorithm import * 

class Rossa(CPAlgorithm):
	"""Rossa's algorithm for finding continuous core-periphery structure.
	
	Examples
	--------
	Create this object.

	>>> import cpalgorithm as cpa	
	>>> rs = cpa.Rossa()
	
	**Core-periphery detection**
	
	Detect core-periphery structure in network G (i.e., NetworkX object):
	
	>>> rs.detect(G) 
	
	Retrieve the ids of the core-periphery pair to which each node belongs:
	
	>>> pair_id = rs.get_pair_id() 
	
	Retrieve the coreness:

	>>> coreness = rs.get_coreness() 
		
	.. note::

	   This algorithm can accept unweighted and weighted networks.
	   The algorithm assigns all nodes into the same core-periphery pair by construction, i.e., c[node_name] =0 for all node_name.
	
	
	.. rubric:: Reference

	F. Rossa, F. Dercole, and C. Piccardi. Profiling core-periphery network structure by random walkers. Scientific Reports, 3, 1467, 2013
	

	"""
	
	
	def __init__(self):
		return
	
	def detect(self, G):
		"""Detect a single core-periphery structure.
	
		Parameters
		----------
		G : NetworkX graph object
		
		Examples
		--------
		>>> import networkx as nx
		>>> import cpalgorithm as cpa
		>>> G = nx.karate_club_graph()  # load the karate club network. 
		>>> rs = cp.Rossa()
		>>> rs.detect(G)
		"""
		
		self.c_, self.x_ = self._detect(G)
		self.Q_ = self._score(G, self.c_, self.x_) 
		self.qs_ = self.Q_
	
	def _score(self, G, c, x):
		nodes = G.nodes()
		xx = np.zeros((len(nodes),1))
		for idx, nd in enumerate(nodes):	
			xx[idx] = x[nd]
		return [1-np.sum(xx) / len(nodes)]

	def _detect(self, G):
		
		node_pairs, w, node2id, id2node = self._to_edge_list(G)

		N = nx.number_of_nodes(G)
		deg = np.array([d[1] for d in G.degree()])
		deg = np.asmatrix(deg)
		M = sum(deg) / 2.0	
		
		A = nx.to_scipy_sparse_matrix(G)
		
		x = np.zeros((N, 1))

		idx = self._argmin2(np.squeeze(np.asarray(deg)))
			
		x[idx] = 1
		ak = deg[0,idx]
		bk = 0
		alpha = np.zeros(N)
		
		for k in range(1, N):
			denom = np.asscalar(np.max([1,np.max(ak * (ak + deg))]) )
			score = (2 * ak * (x.T * A) - bk * deg) / denom
				
			score[ x.T > 0 ] = np.Infinity
			score = np.squeeze(np.asarray(score))
			idx = self._argmin2(score)
			x[idx] = 1
			ak = ak + deg[0, idx]
			bk = np.asscalar(np.dot(x.T* A, x)[0,0])

			alpha[idx] = bk / max(1, ak)
			
		c = dict(zip( [id2node[i] for i in range(N)], np.zeros(N).astype(int)))
		x = dict(zip( [id2node[i] for i in range(N)], alpha))
		return c, x
	
	def _argmin2(self, b): 
		return np.random.choice(np.flatnonzero(b == b.min())) 
