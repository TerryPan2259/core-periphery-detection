/*
*
* Header file of the Divisive algorithm
*
*
* Please do not distribute without contacting the authors.
*
* AUTHOR - Sadamori Kojaku
*
* DATE - 04 July 2018
*/
#ifndef CP_ALGORITHM 
#define CP_ALGORITHM
	#include "cpalgorithm.h" 
#endif

#ifndef BE_ALGORITHM 
#define BE_ALGORITHM
	#include "bealgorithm.h" 
#endif

#include <math.h> 

class Divisive: public CPAlgorithm{
public:
	// Constructor 
	Divisive();
	Divisive(int num_runs);
	
	void detect(const Graph& G);
	
	void calc_Q(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<double>& x,
	    double& Q,
	    vector<double>& q);
	
protected: // function needed to be implemented
	int _num_runs; 
	void _detect_(const Graph& G, vector<int>& c, vector<double>& x, double Q, vector<double>& qs);
	void _louvain(const Graph& G, const int num_of_runs, vector<int>& c, mt19937_64& mtrnd);
	void _louvain_core(const Graph& G, vector<int>& C, const double M, mt19937_64& mtrnd);
	double _calc_dQmod( double d, double degi, double D, double selfw, const double M );
	double _calc_Qmod(const Graph& G, vector<int>& C, const double M);
	void _calc_Q_be(const Graph& G, const vector<int>& c, const vector<double>& x, double& Q, vector<double>& q);
	void _coarsing(const Graph& G, const vector<int>& c,Graph& newG);
	void _modularity_label_switching(const Graph& G,vector<int>& C, const double M,mt19937_64& mtrnd);
	void _subgraph(const Graph& G, vector<bool>&slice, Graph& Gs );	
};



/*-----------------------------
Constructor
-----------------------------*/
Divisive::Divisive(int num_runs):CPAlgorithm(){
	Divisive();
	_num_runs = num_runs;
};

Divisive::Divisive(): CPAlgorithm(){
	_num_runs = 10;
};


/*-----------------------------
Functions inherited from the super class (CPAlgorithm)
-----------------------------*/
void Divisive::calc_Q(
    const Graph& G,
    const vector<int>& c,
    const vector<double>& x,
    double& Q,
    vector<double>& q)
{
    	int K = *max_element(c.begin(), c.end()) + 1;
	vector<double> qtmp(K,0.0);
	q = qtmp;	
	Q=0.0;
	int N = G.get_num_nodes();
	for(int k = 0; k < K; k++){
	
		vector<bool> slice(N,false);
		int Ns = 0;
		for(int i = 0; i < N; i++){
			if(c[i] == k){
				slice[i] = true;
				Ns++;
			}
		}
		int idx = 0;
		vector<double> xs(Ns, 0.0);
		for(int i = 0; i < N; i++){
			if(slice[i]){
				xs[idx] = x[i];
				idx++;
			}
		}
	
		
		
		Graph Gs(0);
		_subgraph(G, slice, Gs);
		
		vector<double> qss;
		double Qs;
		_calc_Q_be(Gs, c, xs, Qs, qss);
		q[k] = Qs;
		Q+=Qs;
	}	
}

void Divisive::detect(const Graph& G){
	_detect_(G, _c, _x, _Q, _q);		
}

            
void Divisive::_detect_(const Graph& G, vector<int>& c, vector<double>& x, double Q, vector<double>& qs){
		
	//Initialise _x and c randomly 
	int N = G.get_num_nodes();
	vector<double> tmp(N, 0.0);
	vector<int> tmpc(N, 0);
	x = tmp;
	c = tmpc;
	uniform_real_distribution<double> dis(0.0, 1.0);	
	
	for(int i = 0;i<N; i++) c[i] = i;
	
	// Detect communities
	_louvain(G, _num_runs, c, _mtrnd); 
	
	// perform the BE algorithm for each community
	BEAlgorithm be = BEAlgorithm(_num_runs);
	
    	int K = *max_element(c.begin(), c.end()) + 1;	
	vector<double> tmp3(K,0.0);
	_q = tmp3;
	_Q = 0;
	for(int k = 0; k < K; k++){
		vector<bool> slice(N,false);
		
		for(int i = 0; i < N; i++){
			if(c[i] == k){
				slice[i] = true;
			}
		}

		Graph Gs(0);
		_subgraph(G, slice, Gs);
		be.detect(Gs);
		
		vector<double> xs = be.get_x();
		vector<double> qss;
		double Qs;
		_calc_Q_be(Gs, c, xs, Qs, qss);
		_q[k] = Qs;
		_Q+=Qs;
				
		int idx = 0;
		for(int i = 0; i < N; i++){
			if(c[i] == k){
				x[i] = xs[idx];
				idx++;
			}
		}
	}
	_c = c;
	_x = x;
} 

void Divisive::_calc_Q_be(
    const Graph& G,
    const vector<int>& c,
    const vector<double>& x,
    double& Q,
    vector<double>& q)
{
	int N = G.get_num_nodes();
	if(N==1){
		Q = 0;
		vector<double> tmp(1,0.0);
		q = tmp;
		return;
	}

	double M = 0.0;
	double pa = 0;
	double pb = 0;
	int nc = 0;
	double mcc = 0;
	for( int i = 0;i < N;i++ ) {
		nc+=(int) x[i];
	
		int sz = G.degree(i);	
		for( int k = 0; k < sz; k++ ) {
			Neighbour nei = G.get_kth_neighbour(i, k);
			int j = nei.get_node(); 
			double w = nei.get_w(); 
			mcc+=w * (x[i]+x[j] - x[i] * x[j]);
			M+=1;
		}
	}
	mcc = mcc/2;
	M = M /2;

	double M_b = (double)(nc * (nc-1) + 2 * nc * (N -nc))/2;
	pa = M / (double)(N * (N-1)/2);
	pb = M_b / (double)(N * (N-1)/2);	
	
	
	Q = ((double)mcc - pa * M_b ) / (sqrt(pa * (1-pa)) * sqrt(pb * (1-pb)) +1e-30);
	Q = Q / (double)(N * (N-1)/2);
		
	if(Q > 1) Q = 1;
	if(Q < -1) Q = -1;
	
	vector<double> qtmp(1,Q);
	q= qtmp;
}

void Divisive::_subgraph(const Graph& G, vector<bool>&slice, Graph& Gs ){
	int N = G.get_num_nodes();
		
	int Ns = 0;
	vector<int> node2id(N,0);
	int idx = 0;
	for(int i = 0; i < N; i ++){
		if(slice[i]){
			node2id[i] = idx;
			Ns++;
			idx+=1;
		}
	}
	
	Graph tt(Ns);
	Gs = tt;
	for(int i =0; i < N;i++){
		int sz = G.degree(i);
		for(int j =0; j < sz;j++){
			Neighbour n = G.get_kth_neighbour(i,j);
			int nei = n.get_node();
			double w = n.get_w();	
			if(slice[i] &  slice[nei]){
				Gs.addEdge(node2id[i], node2id[nei], w);
			}	
		}
	}
}	


double Divisive::_calc_dQmod( double d, double degi, double D, double selfw, const double M ){
	return 2*( d - degi*D/(2.0*M)) + (selfw-degi*degi/(2.0*M)); 
}

double Divisive::_calc_Qmod(
	const Graph& G, 
	vector<int>& C, 
	const double M
	 ){
	
	double retval = 0;
	int N = (int) C.size();
	vector<double> degC(N, 0.0);
	for(int i =0;i<N;i++){
        int di = G.degree(i);
	for(int j =0;j<di;j++){
		Neighbour nei = G.get_kth_neighbour(i,j);
		double w = nei.get_w();
		int k = nei.get_node();

		degC[ C[i] ]+=w;

		if(C[i] == C[k]) {
			retval+=w;	
		}	
	}
	}
	
	for(int i =0;i<N;i++){
		retval-=degC[i]*degC[i]/(2*M);	
	}
	retval/=(2*M);
	return retval;
}




void Divisive::_coarsing(
    	const Graph& G,
    	const vector<int>& c,
    	Graph& newG
	){
		
        int N = (int) c.size();
	
    	int K = *max_element(c.begin(), c.end()) + 1;
	newG = Graph(K);
	for(int i = 0;i<N;i++){
		int mi =  c[i];
		int sz = G.degree(i);
		for(int j = 0;j<sz;j++){
			Neighbour nb = G.get_kth_neighbour(i, j);
			int nei = nb.get_node();
			double w = nb.get_w();
			int mj = c[nei];
			newG.addEdge(mi, mj, w);
		}
	}
	
	newG.compress();
}

void Divisive::_modularity_label_switching(
    	const Graph& G,
	vector<int>& C, 
	const double M,
        mt19937_64& mtrnd
	){
	
        int N= (int)C.size();
	vector<int> ndord(N);
	vector<double> D(N, 0);
	vector<double>deg(N, 0);
	for(int i = 0;i < N;i++) {
		ndord[i] = i;
		deg[i] = G.wdegree(i);
		D[C[i]]+=deg[i];
	};
		
	// --------------------------------
	// Label switching  
	// --------------------------------
	vector<double> toC;
	toC.assign(N,0.0);
	bool isupdated = false;
	int itNum = 0;
	do{
		isupdated = false;
		shuffle(ndord.begin(), ndord.end(),mtrnd);
		for(int i = 0;i <N;i++){
			int nid = ndord[i];
			int ncid = C[nid];
			
			int neighbourNum = G.degree(nid);
			fill(toC.begin(),toC.end(),0.0);	
				
			int newcid = ncid;
			double selfw = 0;
			for(int j = 0;j<neighbourNum;j++){		
				Neighbour nb = G.get_kth_neighbour(nid, j);
				int nei = nb.get_node();
				double w = nb.get_w();

				int cid = C[nei];
	
				if(nid==nei){
					selfw+=w;
				}else{
					toC[cid]+= w;
				}	
			}
			
			double dQold = _calc_dQmod( toC[ncid], deg[nid], D[ncid] - deg[nid], selfw, M );
			double dQ = 0;
			for(int j = 0;j<neighbourNum;j++){
				Neighbour nb = G.get_kth_neighbour(nid, j);
				int nei = nb.get_node();
				int cid = C[nei];
				if(nei==nid) continue;
					
				double dQc=_calc_dQmod( toC[cid], deg[nid], D[cid] - deg[nid]*(double)!!(ncid==cid), selfw, M )-dQold;
			
				if( dQc<dQ ) continue;
					
				newcid = cid;
				dQ = dQc;	
			}
			
			if(dQ< 0) continue;
				
			if( ncid==newcid  ) continue;
			
			
			D[ncid]-=deg[nid];
			D[newcid]+=deg[nid];
			C[nid] = newcid;	
			
			isupdated = true;
		}
		itNum++;
	}while( (isupdated == true) & (itNum<=100) );

	// remove redundant cp
	std::vector<int> labs;
	for(int i=0;i<N;i++){
		int cid = -1;
 		int labsize = (int)labs.size();
		for(int j=0;j<labsize;j++){
			if(labs[j]==C[i]){
				cid = j;
				break;
			}
		}
		
		if (cid<0) {
			labs.push_back(C[i]);
			cid = (int)labs.size()-1;
		}
		C[i] = cid;		
	}
}

void Divisive::_louvain_core(
	const Graph& G, 
	vector<int>& C, 
	const double M,
        mt19937_64& mtrnd
	){
	
	
	
	Graph newG = G; 
	vector<int>Zt = C; 
	vector<int>Ct = C;
	unsigned int prevGraphSize = (int) C.size();
	double Qbest = _calc_Qmod(newG, Zt, M); 

	do{
		prevGraphSize = newG.get_num_nodes();
		
	 	_modularity_label_switching(newG, Zt, M, mtrnd);
		double Qt = _calc_Qmod(newG, Zt, M);
		Graph g;	
		_coarsing(newG, Zt, g);
		newG = g;

		// update C
		// Ct = Ct*Zt;
		int Ctsize = (int) Ct.size();
		for(int i = 0;i<Ctsize;i++){
			Ct[i] = Zt[ Ct[i] ];
		}
		vector<int> tmp(newG.get_num_nodes(),0); 
		for(int i = 0;i<newG.get_num_nodes();i++){
			tmp[i] = i;
		}
		Zt = tmp;
		
		if(Qt>Qbest){
			C = Ct;
			Qbest = Qt;
		}
	}while( newG.get_num_nodes()!= prevGraphSize);
}

void Divisive::_louvain(
    const Graph& G,
    const int num_of_runs,
    vector<int>& C,
    mt19937_64& mtrnd
    ){

    int N = G.get_num_nodes();
    vector<int> tmp(N,0);
    C = tmp;

    double M = 0;
    for(int i =0;i<N;i++){
    	C[i] = i;
	M+=G.wdegree(i);
    }
    M = M / 2; 

    double Q = -1;
    vector<int> cbest;
    for (int i = 0; i < num_of_runs; i++) {
        vector<int> ci = C;
        double Qi = 0.0;

        _louvain_core(G, ci, M, mtrnd);
	Qi = _calc_Qmod(G, ci, M);
	
        if (Qi > Q) {
            Q = Qi;
            cbest = ci;
        }
    }
    C = cbest;
}
