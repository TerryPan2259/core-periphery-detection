import numpy as np
import networkx as nx
import multiprocessing as mp
from scipy.stats import norm

def sz_n(network, c, x):
    return np.bincount(list(c.values())).tolist()

def sz_degree(network, c, x):
    K = max(c.values())+1
    w = [0 for i in range(K)]
    for key, val in c.items():
       w[val]+=network.degree(key)

    return w 

def config_model(G):
    deg = [d[1] for d in G.degree()]
    return nx.configuration_model(deg)	

def erdos_renyi(G):
    n = G.number_of_nodes()
    p = nx.density(G)
    return nx.fast_gnp_random_graph(n, p)

def qstest(c, x, G, cpa, significance_level=0.05, random_net_generator = config_model, sfunc = sz_n, num_of_thread = 4, num_of_rand_net = 500 ):
    q = np.array(cpa.score(G, c, x), dtype = np.float)    
    s = np.array(sfunc(G, c, x) , dtype = np.float)
    C = len(q)

    alpha_corrected = 1.0 - (1.0 - significance_level) ** (1.0 / float(C))
    
    q_tilde = []
    s_tilde = []
    if num_of_thread == 1:
        q_tilde, s_tilde = draw_qs_samples(G, sfunc, cpa, random_net_generator, num_of_rand_net)
    else:
        private_args = [(G, sfunc, cpa, random_net_generator, int(num_of_rand_net / num_of_thread) + 1) for i in range(num_of_thread)]
        pool = mp.Pool(num_of_thread)
        qs_tilde = pool.map(wrapper_draw_qs_samples, private_args)
        for i in range(num_of_thread):
            q_tilde += qs_tilde[i][0] 
            s_tilde += qs_tilde[i][1]
 
    q_tilde = np.array(q_tilde, dtype = np.float)    
    s_tilde = np.array(s_tilde, dtype = np.float)    
    
    q_ave = np.mean(q_tilde)
    s_ave = np.mean(s_tilde)
    q_std = np.std(q_tilde, ddof = 1)
    s_std = np.std(s_tilde, ddof = 1)

    if (s_std <= 1e-30) or (q_std <= 1e-30):
        gamma = 0.0
    else:
        gamma = np.corrcoef(q_tilde, s_tilde)[0, 1]

    h = float(len(q_tilde)) ** (- 1.0 / 6.0)
    p_values = [1.0] * C
    significant = [False] * C
    for cid in range(C):
        if (s_std <= 1e-30) or (q_std <= 1e-30):
            continue    
        w = np.exp(- ( (s[cid] - s_tilde) / (np.sqrt(2.0) * h * s_std) ) ** 2)
        cd = norm.cdf( ( (q[cid] - q_tilde) / (h * q_std) - gamma * (s[cid] - s_tilde) / (h * s_std) ) / np.sqrt(1.0 - gamma * gamma) )    
        denom = sum(w)    
        if denom <= 1e-30:
            continue    
        p_values[cid] = 1.0 - (sum( w * cd ) / denom)
        significant[cid] = p_values[cid] <= alpha_corrected
    return significant, p_values, q_tilde, s_tilde 


# Private function for qstest        
def draw_qs_samples(G, sfunc, cpa, random_net_generator, num_of_rand_net):
    #deg = [x[1] for x in G.degree()]
    q_rand = []
    s_rand = []

    for i in range(num_of_rand_net):
        Gr = random_net_generator(G)
        cpa.detect(Gr)  
        q_rand = q_rand + cpa.score()
        s_rand = s_rand + sfunc(Gr, cpa.get_pair_id(), cpa.is_core()) 
    return q_rand, s_rand


# Private function for qstest        
def wrapper_draw_qs_samples(args):
    return draw_qs_samples(*args)    
