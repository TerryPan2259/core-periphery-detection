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
	"""Rombach's algorithm for finding continuous core-periphery structure.
	
	Parameters
	----------
	num_runs : int
		Number of runs of the algorithm  (optional, default: 1).  

	alpha : float
		Sharpness of core-periphery boundary (optional, default: 0.5). 
	    
	    alpha=0 or alpha=1 gives the fuzziest or sharpest boundary, respectively.   

	beta : float
		Fraction of peripheral nodes (optional, default: 0.8) 

	algorithm : str
		Optimisation algorithm (optional, default: 'ls') 
	    	In the original paper [1], the authors adopted a simulated annealing to optimise the objective function, which is computationally demanding. 
	    	To mitigate the computational cost, a label switching algorithm is implemented in cpalgorithm.
	    	One can choose either algorithm by specifying algorithm='ls' (i.e., label switching) or algorithm='sa' (i.e., simulated annealing).

	.. note::

	   The parameters of the simulated annealing such as the initial temperature and cooling schedule are different from those used in the original paper [1]. 
	
	
	Examples
	--------
	Create this object.

	>>> import cpalgorithm as cpa	
	>>> rb = cpa.Rombach()
	
	**Core-periphery detection**
	
	Detect core-periphery structure in network G (i.e., NetworkX object):
	
	>>> rb.detect(G) 
	
	Retrieve the ids of the core-periphery pair to which each node belongs:
	
	>>> pair_id = rb.get_pair_id() 
	
	Retrieve the coreness:

	>>> coreness = rb.get_coreness() 
		
	.. note::

	   This algorithm can accept unweighted and weighted networks.
	   The algorithm assigns all nodes into the same core-periphery pair by construction, i.e., c[node_name] =0 for all node_name.

	.. rubric:: Reference

        [1] P. Rombach, M. A. Porter, J. H. Fowler, and P. J. Mucha. Core-Periphery Structure in Networks (Revisited). SIAM Review, 59(3):619–646, 2017	

	"""
	
	def __init__(self, num_runs = 10, alpha = 0.5, beta = 0.8, algorithm='ls'):
		self.num_runs = num_runs 
		self.alpha = alpha
		self.beta = beta
		self.algorithm = algorithm
	
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
		>>> rb = cp.Rombach(algorithm='ls') # label switching algorithm
		>>> rb.detect(G)
		>>> rb = cp.Rombach(algorithm='sa') # simulated annealing  
		>>> rb.detect(G)

		"""
		
		Qbest = -100	
		cbest = 0 
		xbest = 0 
		qbest = 0 
		for i in range(self.num_runs):
			if self.algorithm == 'ls':
			
				self._label_switching(G, self.alpha, self.beta)
	
			elif self.algorithm == 'sa':
	
				self._simaneal(G, self.alpha, self.beta)
	
			if Qbest < self.Q_:
				Qbest = self.Q_	
				cbest = self.c_
				xbest = self.x_
				qbest = self.qs_
		
		self.Q_ = Qbest 	
		self.c_ = cbest 
		self.x_ = xbest 
		self.qs_ = qbest 
		
	def _label_switching(self, G, alpha, beta):
		
		node_pairs, w, node2id, id2node = self._to_edge_list(G)

		cppairs = _cp.detect_rombach_ls(edges=node_pairs, ws=w, num_of_runs = 1, alpha = alpha, beta = beta)

		N = len(node2id)	
		self.c_ = dict(zip( [id2node[i] for i in range(N)], cppairs[0].astype(int)))
		self.x_ = dict(zip( [id2node[i] for i in range(N)], cppairs[1]))
		self.Q_ = cppairs[2][0]
		self.qs_ = cppairs[3].tolist()
			
	
	def _simaneal(self, G, alpha, beta):
		
		node_pairs, w, node2id, id2node = self._to_edge_list(G)
		
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
