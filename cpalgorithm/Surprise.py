import _cpalgorithm as _cp
import scipy as sp
import numpy as np
from .CPAlgorithm import * 

class Surprise(CPAlgorithm):
	
	def __init__(self):
		self.num_runs = 1
	
	def detect(self, G):

		nodes = G.nodes()
		N = len(nodes)
		Cbest =np.zeros(N)
		qbest = 0 
		for it in range(self.num_runs):	

			C, q = self._detect(G)
		
			if q < qbest:
				Cbest = C
				qbest = q
		
		# ------------	
		# Post process 
		# ------------
		self.c_ = dict(zip(nodes, np.zeros(N).astype(int)))
		self.x_ = dict(zip(nodes, Cbest.astype(bool)))
		self.Q_ = self._score(G, self.c_, self.x_)[0]
		self.qs_ = [self.Q_]
	
	def _detect(self, G):
		# ----------	
		# Initialise 
		# ----------	
		A = nx.to_scipy_sparse_matrix(nx.Graph(G))
			
		# ------------	
		# Main routine 
		# ------------
		N = A.shape[0]
		C = np.random.randint(2, size=N)
		edges = self._sort_edges(A)
	
		q =  self._calculateSurprise(A, C)

		for itnum in range(5):

			for edge in edges:
				u = edge[0]
				v = edge[1]
			
				# Move node u to the other group and 	
				# compute the changes in the number of edges of different types	
				C[u] = 1-C[u]
			
				# Propose a move
				qp = self._calculateSurprise(A, C)

				if q < qp: # If this is a bad move, then bring back the original
					C[u] = 1-C[u]
				else:
					q = qp
							
						
				if np.random.rand() >0.5:
					# move 1 to 0 for n=3 nodes 
					nnz = np.nonzero(C)[0]
					if nnz.shape[0] ==0:
						continue
					flip = np.random.choice(nnz, np.min([3, nnz.shape[0]]))
					C[flip] = 1 - C[flip]
				else:
					# move 0 to 1 for n=3 nodes 
					nnz = np.nonzero(1-C)[0]
					if nnz.shape[0] ==0:
						continue
					flip = np.random.choice(nnz, np.min([3, nnz.shape[0]]))
					C[flip] = 1 - C[flip]
		
				qp = self._calculateSurprise(A, C)
				if q < qp:
					C[flip] = 1 - C[flip]
				else:
					q = qp

		# Flip group index if group 1 is sparser than group 0
		Nc = np.sum(C)
		Np = N - Nc
		Ecc = np.sum(np.multiply(A @ C, C))
		Epp = np.sum(np.multiply(A @ (1-C), (1-C)))
		if Ecc * (Np *(Np-1)) < Epp * (Nc *(Nc-1)):
			C = 1-C
		return C, q	
	
	def _sort_edges(self, A):
		
		r,c = A.nonzero()
		deg = np.array(A.sum(axis = 0)).ravel()
		
		s = deg[r] * deg[c]
		
		ind = np.argsort(s)
		
		return list(zip(r[ind], c[ind]))

	def _slog(self, numer, denom, s):
		if (s ==0) | (numer < 0) | (denom < 0):
			return 0
	
		denom = denom * (1.0 - 1e-10) + 1e-10
		numer = numer * (1.0 - 1e-10) + 1e-10
	
		v = s * np.log( numer / denom )
		
		return v
	
	def _score(self, G, c, x):

		A = nx.to_scipy_sparse_matrix(nx.Graph(G))
		nodes = G.nodes()
		C = np.array([x[nd] for nd in nodes])
		q = self._calculateSurprise(A, C)
		return [-self._calculateSurprise(A, C)]	
 
	def _calculateSurprise(self, A, x):

		N = A.shape[0]
		Nc = x.sum()
		Np = N - Nc
	
		if (Nc <2) | (Np <2) |(Nc > N-2) | (Np >N-2) :	
			return 0

		L = A.sum()/2
		V = N * (N-1)/2
		Vc = Nc * (Nc-1)/2
		Vcp = Nc * Np 
	
		Ax = A @ x
		lc = np.sum(np.multiply(Ax, x))/2
		lcp = np.sum(np.multiply(Ax, 1-x))
		lp = L - lc - lcp
	
		p = L / (N * (N -1)/2)	
		pc = lc / (Nc * (Nc-1)/2) 
		pp = lp / (Np * (Np-1)/2) 
		pcp = lcp / (Nc * Np)
	
		S = self._slog( p, pp,  2 * L)  +  self._slog((1-p), (1-pp), 2* V-2*L) + self._slog(pp, pc, 2*lc)\
			+ self._slog( (1-pp), (1-pc), 2*Vc  - 2*lc) + self._slog(pp,pcp, 2*lcp)\
			+ self._slog( (1-pp), (1-pcp), 2*Vcp  - 2*lcp)

		return S	
