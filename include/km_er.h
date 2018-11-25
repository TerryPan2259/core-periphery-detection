/*
*
* Header file of the KM-config algorithm (C++ version)
*
*
* An algorithm for finding multiple core-periphery pairs in networks
*
*
* Core-periphery structure requires something else in the network
* Sadamori Kojaku and Naoki Masuda
* Preprint arXiv:1710.07076
* 
*
* Please do not distribute without contacting the authors.
*
*
* AUTHOR - Sadamori Kojaku
*
*
* DATE - 11 Oct 2017
*/
#ifndef CP_ALGORITHM 
#define CP_ALGORITHM
	#include "cpalgorithm.h" 
#endif

class KM_ER: public CPAlgorithm{
public:
	// Constructor 
	KM_ER();
	KM_ER(int num_runs);

	// function needed to be implemented

	void detect(const Graph& G);
	
	void calc_Q(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<double>& x,
	    double& Q,
	    vector<double>& q);
	
protected: 
	
private:
	int _num_runs;
	int _round;
	uniform_real_distribution<double> _udist;
	
	void _km_ER_label_switching(
	    const Graph& G,
	    const int num_of_runs,
	    vector<int>& c,
	    vector<double>& x,
	    double& Q,
	    vector<double>& q,
            mt19937_64& mtrnd
		);

	double _calc_dQ_ER(double d_i_c,
	    double d_i_p,
	    double N_c,
	    double N_p,
	    double selfloop,
	    double x,
	    const double rho);
	
	void _propose_new_label(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<double>& x,
	    const vector<double>& sz_core,
	    const vector<double>& sz_peri,
	    const double rho,
	    const int node_id,
	    int& cprime,
	    double& xprime,
	    double& dQ,
	    mt19937_64& mtrnd
	    );
	
	
	void _km_ER_label_switching_core(
	    const Graph& G,
	    vector<int>& c,
	    vector<double>& x,
	    mt19937_64& mtrnd
	    );

	void _km_ER_louvain(
	    const Graph& G,
	    const int num_of_runs,
	    vector<int>& c,
	    vector<double>& x,
	    double& Q,
	    vector<double>& q,
	    mt19937_64& mtrnd
	    );

/*
	void _km_ER_louvain_core(
		const Graph& G, 
    		vector<int>& c,
    		vector<double>& x,
        	mt19937_64& mtrnd
		);
*/

	void _coarsing(
	    	const Graph& G,
	    	const vector<int>& c,
	    	const vector<double>& x,
	    	Graph& newG,
	    	vector<int>& toLayerId 
		);

	void _relabeling(vector<int>& c);
};


/*-----------------------------
Constructor
-----------------------------*/
KM_ER::KM_ER(int num_runs):CPAlgorithm(){
	KM_ER();
	_num_runs = num_runs;
};

KM_ER::KM_ER():CPAlgorithm(){
	uniform_real_distribution<double> tmp(0.0,1.0);
	_udist = tmp;
	_num_runs = 10;
	_mtrnd = _init_random_number_generator();
};


/*-----------------------------
Functions inherited from the super class (CPAlgorithm)
-----------------------------*/
void KM_ER::detect(const Graph& G){
	//_km_ER_label_switching(G, _num_runs, _c, _x, _Q, _q, _mtrnd);
	_km_ER_louvain(G, _num_runs, _c, _x, _Q, _q, _mtrnd);
	//_km_ER_label_switching(G, _num_runs, _c, _x, _Q, _q, _mtrnd);
}

void KM_ER::calc_Q(
    const Graph& G,
    const vector<int>& c,
    const vector<double>& x,
    double& Q,
    vector<double>& q)
{
    int N = G.get_num_nodes();
    int K = *max_element(c.begin(), c.end()) + 1;
    q.assign(K, 0.0);
    vector<double> Nc(K, 0.0);
    vector<double> Np(K, 0.0);

    double double_M = 0.0;
    for (int i = 0; i < N; i++) {
	int sz = G.degree(i);
        for (int j = 0; j < sz; j++) {
	    Neighbour nn = G.get_kth_neighbour(i, j);
	    int nei = nn.get_node();
	    double wj = nn.get_w();
            q[c[i]] += wj * !!(c[i] == c[nei]) * (x[i] + x[nei] - x[i] * x[nei]);
	    double_M += wj;
        }
	Nc[c[i]]+=x[i];
	Np[c[i]]+=(1-x[i]);
    }
    Q = 0;
    double rho = double_M / (double)(N * N); 
    for (int k = 0; k < K; k++) {
        q[k] = (q[k] - rho * (Nc[k]*Nc[k] + 2 * Nc[k] * Np[k] )) / double_M;
        Q += q[k];
    }
}

/*-----------------------------
Private functions (internal use only)
-----------------------------*/
void KM_ER::_km_ER_label_switching(
    const Graph& G,
    const int num_of_runs,
    vector<int>& c,
    vector<double>& x,
    double& Q,
    vector<double>& q,
    mt19937_64& mtrnd
    )
{
    int N = G.get_num_nodes();
    c.clear();
    x.clear();
    c.assign(N, 0);
    x.assign(N, 1);

    Q = -1;
    for (int i = 0; i < num_of_runs; i++) {
        vector<int> ci;
        vector<double> xi;
        vector<double> qi;
        double Qi = 0.0;

        _km_ER_label_switching_core(G, ci, xi, _mtrnd);

        calc_Q(G, ci, xi, Qi, qi);
        if (Qi > Q) {
	    c = ci;
	    x = xi;
	    q.clear();
	    q = qi;
            Q = Qi;
        }
    }
}
/*
void KM_ER::_km_ER_label_switching(
    const Graph& G,
    const int num_of_runs,
    vector<int>& c,
    vector<double>& x,
    double& Q,
    vector<double>& q,
    mt19937_64& mtrnd
    )
{
    int N = G.get_num_nodes();
    c.clear();
    x.clear();
    c.assign(N, 0);
    x.assign(N, true);

    // Generate \hat q^{(s)} and \hat n^{(s)} (1 \leq s \leq S) 
    // create random number generator per each thread
    int numthread = 1;
    #ifdef _OPENMP
    	# pragma omp parallel
    	{
    		numthread = omp_get_num_threads();
    	}
    #endif
    cout<<numthread<<endl;
    vector<mt19937_64> mtrnd_list(numthread);
    for(int i = 0; i < numthread; i++){
	mt19937_64 mtrnd = _init_random_number_generator();
	mtrnd_list[i] = mtrnd;
    }

    Q = -1;
    #ifdef _OPENMP
    #pragma omp parallel for shared(c, x, Q, q, N, G, mtrnd_list)
    #endif
    for (int i = 0; i < num_of_runs; i++) {
        vector<int> ci;
        vector<double> xi;
        vector<double> qi;
        double Qi = 0.0;

        int tid = 0;
    	#ifdef _OPENMP
        	tid = omp_get_thread_num();
    	#endif
	
        mt19937_64 mtrnd = mtrnd_list[tid];
        _km_ER_label_switching_core(G, ci, xi, mtrnd);

        calc_Q(G, ci, xi, Qi, qi);

        #pragma omp critical
        {
        	if (Qi > Q) {
		    c = ci;
		    x = xi;
		    q.clear();
		    q = qi;
        	    Q = Qi;
        	}
	}
    }
}
*/

double KM_ER::_calc_dQ_ER(double d_i_c,
    double d_i_p,
    double N_c,
    double N_p,
    double selfloop,
    double x,
    const double rho)
{
    return 2 * (d_i_c + d_i_p * (x) - rho * (N_c  + N_p * x) ) + x * (selfloop - 2 * rho);
}


void KM_ER::_propose_new_label(
    const Graph& G,
    const vector<int>& c,
    const vector<double>& x,
    const vector<double>& sz_core,
    const vector<double>& sz_peri,
    const double rho,
    const int node_id,
    int& cprime,
    double& xprime,
    double& dQ,
    mt19937_64& mtrnd)
{
    int N = G.get_num_nodes();
    int neighbourNum = G.degree(node_id);

    vector<double> edges_to_core(N, 0.0);
    vector<double> edges_to_peri(N, 0.0);

    double selfloop = 0;
    for (int j = 0; j < neighbourNum; j++) {
	Neighbour nn = G.get_kth_neighbour(node_id, j);
	int nei = nn.get_node();
	double wj = nn.get_w();
	
	if(node_id == nei){
		selfloop+= wj;
		continue;
	}
	
        edges_to_core[c[nei]] += wj * x[nei];
        edges_to_peri[c[nei]] += wj * (1-x[nei]);
    }

    double N_core = sz_core[c[node_id]] -  x[node_id];
    double N_peri = sz_peri[c[node_id]] -  (1-x[node_id]);

    double dQold = _calc_dQ_ER(edges_to_core[c[node_id]], edges_to_peri[c[node_id]],
        N_core, N_peri, selfloop, x[node_id], rho);

    dQ = 0;
    for (int j = 0; j < neighbourNum; j++) {
	Neighbour nn = G.get_kth_neighbour(node_id, j);
	int nei = nn.get_node();
	//double wj = nn.get_w();

        int cid = c[nei];

        N_core = sz_core[cid] - x[node_id] * (double)!!( c[node_id] == cid);
        N_peri = sz_peri[cid] - (1-x[node_id]) * (double)!!(c[node_id] == cid);

        double Q_i_core = _calc_dQ_ER(edges_to_core[cid], edges_to_peri[cid],
            N_core, N_peri, selfloop, 1, rho);
        double Q_i_peri = _calc_dQ_ER(edges_to_core[cid], edges_to_peri[cid],
            N_core, N_peri, selfloop, 0, rho);
        Q_i_core -= dQold;
        Q_i_peri -= dQold;

        if (MAX(Q_i_core, Q_i_peri) < dQ)
            continue;

        if (Q_i_peri < Q_i_core) {
            xprime = 1.0;
            cprime = cid;
            dQ = Q_i_core;
        }
        else if (Q_i_peri > Q_i_core) {
            xprime = 0.0;
            cprime = cid;
            dQ = Q_i_peri;
        }
        else {
            cprime = cid;
	    if(_udist(mtrnd) < 0.5){
		xprime = 1;
	    }else{
		xprime = 0;
	    }
            dQ = Q_i_core;
        }
    }
}

void KM_ER::_km_ER_label_switching_core(
    const Graph& G,
    vector<int>& c,
    vector<double>& x,
    mt19937_64& mtrnd
    )
{
    /* Variable declarations */
    int N = G.get_num_nodes();
    vector<double> sz_core(N, 0.0);
    vector<double> sz_peri(N, 0.0);
    vector<int> order(N);
    double M = 0;
    bool isupdated = false;
    c.clear();
    x.clear();
    c.assign(N, 0);
    x.assign(N, 1.0);
    _round = 0;
    for (int i = 0; i < N; i++) {
        order[i] = i;
        c[i] = i;
        sz_core[i] += x[i];
        sz_peri[i] += 1-x[i];
        M += G.wdegree(i);
    };
    _round++;
    double rho = M / (double)( N * N );
    M = M / 2;

    /* Label switching algorithm */
    do {
        isupdated = false;
        shuffle(order.begin(), order.end(), mtrnd);

        for (int scan_count = 0; scan_count < N; scan_count++) {
            int i = order[scan_count];

            int cprime = c[i]; // c'
            double xprime = x[i]; // x'

            double dQ = 0;
            _propose_new_label(G, c, x, sz_core, sz_peri,
                rho, i, cprime, xprime, dQ, mtrnd);

            if (dQ <= 0)
                continue;
		
            if ( (c[i] == cprime) & (x[i] == xprime) )
                continue;


            sz_core[c[i]] -= x[i];
            sz_peri[c[i]] -= 1-x[i];

            sz_core[cprime] += xprime;
            sz_peri[cprime] += (1-xprime);

            c[i] = cprime;
            x[i] = xprime;

            isupdated = true;
        }
	_round++;
    } while (isupdated == true);

    /* Remove empty core-periphery pairs */
    std::vector<int> labs;
    for (int i = 0; i < N; i++) {
        int cid = -1;
	int labsize = (int) labs.size();
        for (int j = 0; j < labsize; j++) {
            if (labs[j] == c[i]) {
                cid = j;
                break;
            }
        }

        if (cid < 0) {
            labs.push_back(c[i]);
            cid = (int) labs.size() - 1;
        }
        c[i] = cid;
    }
}

/* Louvain algorithm */
void KM_ER::_km_ER_louvain(
	const Graph& G, 
	const int num_of_runs,
    	vector<int>& c,
    	vector<double>& x,
	double& Q,
	vector<double>& q,
        mt19937_64& mtrnd
	){

	// Intiialise variables	
	int N = G.get_num_nodes();
	c.clear();
	x.clear();
    	c.assign(N, 0); 
    	x.assign(N, 1);
    	for (int i = 0; i < N; i++) c[i] = i;
	
	vector<int>ct = c; // label of each node at tth iteration
	vector<double>xt = x; // label of each node at tth iteration. 
	Graph cnet_G; // coarse network
	vector<int> toLayerId; //toLayerId[i] maps 2*c[i] + x[i] to the id of node in the coarse network 
	_coarsing(G, ct, xt, cnet_G, toLayerId); // Initialise toLayerId

	Q = 0; // quality of the current partition

	int cnet_N;
	do{
		cnet_N = cnet_G.get_num_nodes();
		
		// Core-periphery detection	
		vector<int> cnet_c; // label of node in the coarse network, Mt 
		vector<double> cnet_x; // label of node in the coarse network, Mt 
		double Qt = 0; vector<double> qt;
		_km_ER_label_switching(cnet_G, num_of_runs, cnet_c, cnet_x, Qt, qt, mtrnd);
		//_km_ER_label_switching_core(cnet_G, cnet_c, cnet_x, mtrnd);
	
		// Update the label of node in the original network, ct and xt.	
		for(int i = 0; i< N; i++){
			int cnet_id = toLayerId[2 * ct[i] + (int)xt[i]];
			ct[i] = cnet_c[ cnet_id ];
			xt[i] = cnet_x[ cnet_id ];
		}
 		
		// Compute the quality       	
		//calc_Q(G, ct, xt, Qt, qt);
		calc_Q(cnet_G, cnet_c, cnet_x, Qt, qt);

		if(Qt>=Q){ // if the quality is the highest among those detected so far
			c = ct;
			x = xt;
			Q = Qt;
			q = qt;
		}
	
		// Coarsing	
		Graph new_cnet_G; 
		_coarsing(cnet_G, cnet_c, cnet_x, new_cnet_G, toLayerId);
		cnet_G = new_cnet_G;
		
		
		int sz = cnet_G.get_num_nodes();
		if(sz == cnet_N) break;	
			
	}while( true );

	_relabeling(c);
}

void KM_ER::_coarsing(
    	const Graph& G,
    	const vector<int>& c,
    	const vector<double>& x,
    	Graph& newG,
    	vector<int>& toLayerId 
	){
		
        int N = (int) c.size();
	vector<int> ids(N,0);
    	int maxid = 0;
	for(int i = 0;i<N;i++){
		ids[i] = 2 * c[i] + (int)x[i];
		maxid = MAX(maxid, ids[i]);
	}
	_relabeling(ids);
	toLayerId.clear();
	toLayerId.assign(maxid+1,0);
	for(int i = 0;i<N;i++){
		toLayerId[2 * c[i] + (int)x[i]] = ids[i];
	}
	
	
    	int K = *max_element(ids.begin(), ids.end()) + 1;
	newG = Graph(K);
	for(int i = 0;i<N;i++){
		int mi = 2 * c[i] + (int)x[i];
		int sz = G.degree(i);
		for(int j = 0;j<sz;j++){
			Neighbour nb = G.get_kth_neighbour(i, j);
			int nei = nb.get_node();
			double w = nb.get_w();
				
			int mj = 2 * c[nei] + (int)x[nei];

			int sid = toLayerId[mi];
			int did = toLayerId[mj];
			newG.addEdge(sid, did, w);
		}
	}
	
	newG.compress();
}

void KM_ER::_relabeling(
    	vector<int>& c
	){

    int N = (int)c.size(); 
    std::vector<int> labs;
    for (int i = 0; i < N; i++) {
        int cid = -1;
	int labsize = (int)labs.size();
        for (int j = 0; j < labsize; j++) {
            if (labs[j] == c[i]) {
                cid = j;
                break;
            }
        }

        if (cid < 0) {
            labs.push_back(c[i]);
            cid = (int) labs.size() - 1;
        }
        c[i] = cid;
    }
}
