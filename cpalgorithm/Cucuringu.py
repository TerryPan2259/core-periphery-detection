import _cpalgorithm as _cp
from .CPAlgorithm import * 

class Cucuringu_low_rank_core(CPAlgorithm):
	
	def __init__(self):
		self.num_runs = 10 
		sefl.beta = 0.1
	
	def detect(self, G):

		self.c_ = dict(zip( [id2node[i] for i in range(N)], cppairs[0].astype(int)))
		self.x_ = dict(zip( [id2node[i] for i in range(N)], cppairs[1].astype(int)))
		self.Q_ = cppairs[2][0]
		self.qs_ = cppairs[3].tolist()
	
	def _score(self, G, c, x):

		return result[1].tolist()	
	
	def significance(self):
		return self.pvalues
	
	def find_cut(G, score, b):

		N = nx.number_of_nodes(G)
		M = nx.number_of_edges(G) 
		A = nx.to_scipy_sparse_matrix(G)
		q = np.zeros(N)
		deg = np.array([d[1] for d in G.degree()])
		od = (-deg).argsort()

		for i in range(b, N-b):
			x = np.zeros((N,1))
			x[od[0:i]] = 1
			
			Mcc = np.dot(x.T * A, x)/2
			Mcp = np.dot(x.T * A, (1-x))
			Mpp = np.dot(x.T * A, x)/2
			
			q[i] = Mcc/(i*(i-1)/2) + Mcp/(i*(N-i)) - Mpp / ((N-i)*((N-i)-1)/2)
		
		idx = np.argmax(q)
		Q = q[idx]
		Q = Q/N
		
		x = np.zeros((N,1))
		x[od[0:idx]] = 1
		c = dict(zip( [id2node[i] for i in range(N)], np.zeros(N)))
		x = dict(zip( [id2node[i] for i in range(N)], x.astype(int)))
		
		return c, x
	
	def lapsgn_core(G):
		N = nx.number_of_nodes(G)
		M = nx.number_of_edges(G) 
		A = nx.to_scipy_sparse_matrix(G)
		deg = np.array([d[1] for d in G.degree()])
		denom =  1.0 / deg
		denom[np.isnan(denom)] = 0
		T = diags(denom) * A
		
			
		v, d = eigs(T, k=1, sigma=1e-2);
		
		v = sign(v);
		c = dict(zip( [id2node[i] for i in range(N)], np.zeros(N)))
		x = dict(zip( [id2node[i] for i in range(N)], int(np.sign(v)>0)))
		return c, x	
