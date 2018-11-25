import _cpalgorithm as _cp
import scipy as sp
import numpy as np
from .CPAlgorithm import * 

class Surprise(CPAlgorithm):
	
	def __init__(self):
		self.num_runs = 10
	
	def detect(self, G):
		
		C = self._detect(G)
		# ------------	
		# Post process 
		# ------------
		self.c_ = dict(zip( [id2node[i] for i in range(N)], cppairs[0].astype(int)))
		self.x_ = dict(zip( [id2node[i] for i in range(N)], cppairs[1].astype(bool)))
		self.Q_ = cppairs[2][0]
		self.qs_ = cppairs[3].tolist()
	
	def _detect(self, G):
		# ----------	
		# Initialise 
		# ----------	
		A = nx.to_scipy_sparse_matrix(G)
			
		
		# ------------	
		# Main routine 
		# ------------
		N = A.shape[0]
		C = np.random.randint(2, size=(N, 1))
		edges = self._sort_edges(A)
		q = self._calculateSurprise(A, C)

		for itnum in range(10):

			for edge in edges:
				u = edge[0]
				v = edge[1]
	
				C[u] = 1-C[u]
				
				# Propose a move
				qp = self._calculateSurprise(A, C)

				if q < qp: # If this is a bad move, then bring back the original
					C[u] = 1-C[u]
				else:
					q = qp
							
					
				if np.random.rand() >0:
					# move 1 to 0 for n=3 nodes 
					nnz = np.nonzero(C)[0]
					flip = np.random.choice(nnz, np.min([3, nnz.shape[0]]))
					C[flip] = 1 - C[flip]
				else:
					# move 0 to 1 for n=3 nodes 
					nnz = np.nonzero(1-C)[0]
					flip = np.random.choice(nnz, np.min([3, nnz.shape[0]]))
					C[flip] = 1 - C[flip]
		
				qp = self._calculateSurprise(A, C)
				if q < qp:
					C[flip] = 1 - C[flip]
				else:
					q = qp
		return C	

	def _sort_edges(self, A):
		
		r,c = A.nonzero()
		deg = np.array(A.sum(axis = 0)).ravel()
		
		s = deg[r] * deg[c]
		
		ind = np.argsort(s)
		
		return list(zip(r[ind], c[ind]))
		
	
		
	def _calculateAndUpdateSurprise(self, A, C, Cp):

		sp = self._calculateSurprise(A, Cp)
		if sp < s:
			return Cp, sp 
		else:
			return C, s
		
		
	def _calculateSurprise(self, A, C):
		N = A.shape[0]
		Nc = C.sum()
		Np = N - Nc
	
		if (Nc == 0) | (Np ==0) |(Nc == N) | (Np ==N) :	
			return 1e+30

		L = A.sum()/2
		V = N * (N-1)/2
		Vc = Nc * (Nc-1)/2
		Vcp = Nc * Np 
		
		lc = (C.T @ A @ C)[0,0]/2
		lcp = ((1-C.T) @ A @ C)[0,0]
		lp = L - lc - lcp
	
	
		p = L / (N * (N -1)/2)	
		pc = lc / (Nc * (Nc-1)) 
		pp = lp / (Np * (Np-1)) 
		pcp = lcp / (Nc * Np)
	
#		S = np.power( p / pp, L) * np.power( (1-p)/(1-pp), V-L) * np.power(pp/pc, lc)\
#			* np.power( (1-pp)/ (1-pc), Vc  - lc) * np.power(pp/pcp, lcp)\
#			* np.power( (1-pp)/ (1-pcp), Vcp  - lcp)
		
		S = self._slog( p, pp,  2 * L)  +  self._slog((1-p), (1-pp), 2* V-2*L) + self._slog(pp, pc, 2*lc)\
			+ self._slog( (1-pp), (1-pc), 2*Vc  - 2*lc) + self._slog(pp,pcp, 2*lcp)\
			+ self._slog( (1-pp), (1-pcp), 2*Vcp  - 2*lcp)
		return S	

	def _slog(self, numer, denom, s):
		if s ==0:
			return 0
		v = s * np.log( numer / (denom+1e-30) )
		
		return v
	
	def _score(self, G, c, x):

		return 

	def significance(self):
		return self.pvalues	
