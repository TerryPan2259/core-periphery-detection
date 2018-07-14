import _cpalgorithm as _cp
from simanneal import Annealer
import random
from .CPAlgorithm import * 

class SimAlg(Annealer):
	def __init__(self, A, x, alpha, beta):
		self.state=x
		self.A = A
		self.alpha = alpha
		self.beta = beta
		self.Tmax = 1
		self.Tmin = 1e-8 
		self.steps = 10000 
		self.updates = 100	

	def move(self):
		"""Swaps two nodes"""
		a = random.randint(0, len(self.state) - 1)
		b = random.randint(0, len(self.state) - 1)
		self.state[[a,b]] = self.state[[b,a]]
		
	def energy(self):
		return self.eval(self.state)		
	
	def eval(self, od):
		
		x = self.corevector(od, self.alpha, self.beta)
		
		return -np.asscalar(np.dot(x.T * self.A, x)[0,0])
	

	def corevector(self, x, alpha, beta):
		N = len(x);
		bn = np.floor(beta * N)
		cx = (x<=bn).astype(int)
		px = (x>bn).astype(int)

		c = (1.0 - alpha) / (2.0 * bn) * x * cx + ((x * px - bn ) * (1.0 - alpha) / (2.0 * (N - bn)) + (1.0 + alpha) / 2.0) * px;
		return np.asmatrix(np.reshape(c, (N,1)))

class Rombach(CPAlgorithm):
	
	def __init__(self):
		self.num_runs = 10 
		self.alpha = 0.5
		self.beta = 0.8
		self.algorithm = 'ls'
	
	def detect(self, G, alpha = -1, beta = -1, algorithm = []):
		if alpha >0:	
			self.alpha = alpha		

		if beta >0:	
			self.beta = beta	
		
		if beta >0:	
			self.beta = beta	
		
		if algorithm:
			self.algorithm = algorithm
		
		if self.algorithm == 'ls':	
			self.label_switching(G, self.alpha, self.beta)
		elif self.algorithm == 'sa':
			self.simaneal(G, self.alpha, self.beta)
		
		
	def label_switching(self, G, alpha, beta):
		
		node_pairs, w, node2id, id2node = self.to_edge_list(G)

		cppairs = _cp.detect_rombach_ls(edges=node_pairs, ws=w, num_of_runs = self.num_runs, alpha = alpha, beta = beta)

		N = len(node2id)	
		self.c_ = dict(zip( [id2node[i] for i in range(N)], cppairs[0].astype(int)))
		self.x_ = dict(zip( [id2node[i] for i in range(N)], cppairs[1]))
		self.Q_ = cppairs[2][0]
		self.qs_ = cppairs[3].tolist()
			
	
	def simaneal(self, G, alpha, beta):
		
		node_pairs, w, node2id, id2node = self.to_edge_list(G)
		
		A = nx.to_scipy_sparse_matrix(G)
		N = nx.number_of_nodes(G)

		nodes = list(range(N))
		random.shuffle(nodes)
		nodes = np.array(nodes)
		
		sa = SimAlg(A, nodes, self.alpha, self.beta)
	
		od, self.Q_ = sa.anneal()
		self.Q_ *= -1
		
		x = sa.corevector(od, self.alpha, self.beta)
		x = x.T.tolist()[0]
		
		self.c_ = dict(zip( [id2node[i] for i in range(N)], np.zeros(N)))
		self.x_ = dict(zip( [id2node[i] for i in range(N)], x))
		self.qs_ = [-sa.eval(od)]
	
	def _score(self, G, c, x):
		A = nx.to_scipy_sparse_matrix(G)
				
		N = nx.number_of_nodes(G)
		xx = np.zeros((N,1))
		for idx, nd in enumerate(G.nodes()):	
			xx[idx] = x[nd]
		
		return [np.asscalar(np.dot(xx.T * A, xx)[0,0])]
	
	def significance(self):
		return self.pvalues
	

