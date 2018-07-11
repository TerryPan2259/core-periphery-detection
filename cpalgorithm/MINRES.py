import _cpalgorithm as _cp
from .CPAlgorithm import * 

class MINRES(CPAlgorithm):
	
	def __init__(self):
		self.num_runs = 0 
	
	def detect(self, G):

		node_pairs, w, node2id, id2node = self.to_edge_list(G)

		cppairs = _cp.detect_minres(edges=node_pairs, ws=w)
		
		N = len(id2node) 
		self.c_ = dict(zip( [id2node[i] for i in range(N)], cppairs[0].astype(int)))
		self.x_ = dict(zip( [id2node[i] for i in range(N)], cppairs[1].astype(bool)))
		self.Q_ = cppairs[2][0]
		self.qs_ = cppairs[3].tolist()
	
	def _score(self, G, c, x):

		node_pairs, w, node2id, id2node = self.to_edge_list(G)
	
		N = len(id2node)
		_c = np.array([ c[id2node[i]]  for i in range(N) ])
		_x = np.array([ x[id2node[i]]  for i in range(N) ])
	
		result = _cp.calc_Q_minres(edges=node_pairs, ws=w, c=_c, x=_x)
		return result[1].tolist()
	
	def significance(self):
		return self.pvalues	
