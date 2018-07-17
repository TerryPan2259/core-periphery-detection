"""Detect core-periphery structure in empirical networks.

"""
import cpalgorithm as cp
import scipy.stats as stats
import network as nx

def profiling_cp_structure(G, alg):
	
	# cp detection
	alg.detect(G) 

	# retrieve result	
	pair_id = alg.get_core_id()
	coreness = alg.get_coreness()
	
	# statistical test	
	pair_id, coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, alg)
	
	# Number of pairs
	C = np.max(pair_id) + 1
	
	# Average coreness value (in case of disc. cp, this is the fraction of core nodes)
	fc = np.mean(coreness)
	
	#
		
		
	# compute the correlation between degree and the coreness 
	deg = np.array([d[1] for d in G.degree()])
	corr = stats.kendalltau(deg, coreness)
	
	
	
		
	
	
		
	
	
