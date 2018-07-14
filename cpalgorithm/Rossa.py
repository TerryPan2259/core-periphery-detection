import _cpalgorithm as _cp
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from .CPAlgorithm import * 

class Rossa(CPAlgorithm):
	
	def __init__(self):
		self.num_runs = 10 
		self.beta = 0.1
	
	def detect(self, G):
		
		self.c_, self.x_ = self._detect(G)
		self.Q_ = self._score(G, self.c_, self.x_) 
		self.qs_ = self.Q_
	
	def _score(self, G, c, x):
		nodes = G.nodes()
		xx = np.zeros((len(nodes),1))
		for idx, nd in enumerate(nodes):	
			xx[idx] = x[nd]
		return [2 * (np.sum(xx) - np.max(xx)) / max(1, len(xx) -2)]

	def significance(self):
		return self.pvalues

	def _detect(self, G):
		
		node_pairs, w, node2id, id2node = self.to_edge_list(G)

		N = nx.number_of_nodes(G)
		deg = np.array([d[1] for d in G.degree()])
		deg = np.asmatrix(deg)
		M = sum(deg) / 2.0	
		
		A = nx.to_scipy_sparse_matrix(G)
		
		x = np.zeros((N, 1))

		idx = self.argmin2(np.squeeze(np.asarray(deg)))
			
		x[idx] = 1
		ak = deg[0,idx]
		bk = 0
		alpha = np.zeros(N)
		
		for k in range(1, N):
			denom = np.asscalar(np.max([1,np.max(ak * (ak + deg))]) )
			score = (2 * ak * (x.T * A) - bk * deg) / denom
				
			score[ x.T > 0 ] = np.Infinity
			score = np.squeeze(np.asarray(score))
			idx = self.argmin2(score)
			x[idx] = 1
			ak = ak + deg[0, idx]
			bk = np.asscalar(np.dot(x.T* A, x)[0,0])

			alpha[idx] = bk / max(1, ak)
			
		c = dict(zip( [id2node[i] for i in range(N)], np.zeros(N)))
		x = dict(zip( [id2node[i] for i in range(N)], alpha))
		return c, x
	
	def argmax2(self, b): 
		return np.random.choice(np.flatnonzero(b == b.max())) 
	
	def argmin2(self, b): 
		return np.random.choice(np.flatnonzero(b == b.min())) 
