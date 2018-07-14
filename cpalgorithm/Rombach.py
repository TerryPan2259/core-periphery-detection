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

	def move(self):
		"""Swaps two nodes"""
		a = random.randint(0, len(self.state) - 1)
		b = random.randint(0, len(self.state) - 1)
		self.state[[a,b]] = self.state[[b,a]]
		
	def energy(self):
		
		x = self.corevector(self.state, self.alpha, self.beta)
		
		return -np.dot(x.T * self.A, x)[0,0]
	

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
	
	def detect(self, G):
		
		A = nx.to_scipy_sparse_matrix(G)
		N = nx.number_of_nodes(G)

		nodes = list(range(N))
		random.shuffle(nodes)
		nodes = np.array(nodes)
		
		sa = SimAlg(A, nodes, self.alpha, self.beta)
	
		od, self.Q_ = sa.anneal()
		
		x = sa.core_vector(od, self.alpha, self.beta)
		
		c = dict(zip( [id2node[i] for i in range(N)], np.zeros(N)))
		print(x.flatten) 
		self.x_ = dict(zip( [id2node[i] for i in range(N)], x.flatten()))
		print(od)	
			
			
	
	def _score(self, G, c, x):
		return 

	
	def significance(self):
		return self.pvalues
	

