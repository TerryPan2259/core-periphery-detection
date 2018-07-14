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

class Rombach_LS: public CPAlgorithm{
public:
	// Constructor 
	Rombach_LS();
	Rombach_LS(int num_runs);
	Rombach_LS(int num_runs, double alpha, double beta);

	// function needed to be implemented
	void detect(const Graph& G);
	
	void calc_Q(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<double>& x,
	    double& Q,
	    vector<double>& q);

protected:
	int _num_runs;
	int _N;
	double _alpha;
	double _beta;

	double _calc_coreness(int ord);
	double _rowSum_score(const Graph& G, vector<int>& ndord, int nid, int sid);
	double _calc_swap_gain(const Graph& G, vector<int>& ndord, int nid, int sid );
	void _swap_nodes(vector<int>&ndord, int nid, int nextnid );
	void _rombach_label_switching(const Graph& G, vector<int>& ords);
};


/*-----------------------------
Constructor
-----------------------------*/
Rombach_LS::Rombach_LS():CPAlgorithm(){
};
Rombach_LS::Rombach_LS(int num_runs):CPAlgorithm(){
	_num_runs = num_runs;
};
Rombach_LS::Rombach_LS(int num_runs, double alpha, double beta):CPAlgorithm(){
	_num_runs = num_runs;
	_alpha = alpha;
	_beta = beta;
};

/*---------------------------
Functions inherited from the super class (CPAlgorithm)
-----------------------------*/
void Rombach_LS::detect(const Graph& G){

	vector<int> ords;
	_N = G.get_num_nodes();

	_rombach_label_switching(G, ords);

	vector<int> tmp2(_N, 0);	
	_c = tmp2;
	
	vector<double> tmp(_N, 0.0);	
	_x = tmp;	
	for(int i = 0; i < _N; i++){
		_x[i] = _calc_coreness(ords[i]);
	}	
}

void Rombach_LS::calc_Q(
    const Graph& G,
    const vector<int>& c,
    const vector<double>& x,
    double& Q,
    vector<double>& q)
{

	Q = 0;
	int N = G.get_num_nodes();
	for(int i = 0; i < N; i++){
		int sz = G.degree(i);
        	for (int j = 0; j < sz; j++) {
			Neighbour nn = G.get_kth_neighbour(i, j);
			int nei = nn.get_node();
		    	double wj = nn.get_w();
			Q+= wj * x[i] * x[nei];
        	}
	}
	
	vector<double> tmp(1,0);
	tmp[0] = Q;
	q = tmp;
}


/*-----------------------------
Private functions (internal use only)
-----------------------------*/
double Rombach_LS::_calc_coreness( int ord){
	double c = 0.0;
	double bn = floor( _N * _beta);
	if(ord <= bn ){
		c = (1.0-_alpha) / (2.0*bn) * (double)ord;
	}else{
		c = ((double)ord-bn)*(1.0-_alpha)/(2*((double)_N-bn))  + (1+_alpha)/2.0;
	}	
	return c;
}


double Rombach_LS::_rowSum_score(const Graph& G, vector<int>&ndord, int nid, int sid){
	double retval = 0;
	int sz = G.degree(nid);
	for(int k = 0; k < sz; k++){
		Neighbour n = G.get_kth_neighbour(nid, k);
		int id = n.get_node();
		double w = n.get_w();
		//for( auto id : adjList[nid] ){
		if(id==sid) continue;
		retval+=w * _calc_coreness( ndord[id] );
	}
	return retval;
}


double Rombach_LS::_calc_swap_gain(const Graph& G,  vector<int>&ndord, int nid, int sid ){
	double c_nid = _calc_coreness( ndord[nid] );
	double c_sid = _calc_coreness( ndord[sid] );
	double rowSum_nid = _rowSum_score(G, ndord, nid, sid);
	double rowSum_sid = _rowSum_score(G, ndord, sid, nid);
	double dQ = (c_sid - c_nid) * rowSum_nid + (c_nid - c_sid) * rowSum_sid;
	return dQ; 
}
void Rombach_LS::_swap_nodes(vector<int>&ndord, int nid, int nextnid ){
	// update ndord;	
	int tmp = ndord[nid];
	ndord[nid] = ndord[nextnid];
	ndord[nextnid]  = tmp;
	
}

void Rombach_LS::_rombach_label_switching(const Graph& G, vector<int>& nds){

	// Initialise	
	int N = G.get_num_nodes();
	vector<int> ndord(N,0);
	for(int i = 0;i < N;i++) ndord[i] = i;
	vector<int> ord = ndord;
	shuffle(ndord.begin(), ndord.end(), _mtrnd);

	// Label Switching algorithm	
	bool isupdated = false;
	int itmax = 100;int itnum=0;	
	do{
		isupdated = false;	
		shuffle(ord.begin(), ord.end(), _mtrnd);
		for(int i = 0; i < _N; i++){	
			int nid = ord[i]; // Get the id of node we shall update

			// calculate gain obtained by swapping label of nid and other labls
			int nextnid = nid; 

			double dQmax = -_N;

			for( int sid =0;sid<_N;sid++ ){
				if(sid==nid) continue;
				double dQ = _calc_swap_gain(G, ndord, nid, sid);
				if(dQmax < dQ){
					nextnid = sid;
					dQmax = dQ;
				}
			}
			if( dQmax <= std::numeric_limits<double>::epsilon()) continue;
			isupdated = true;	
			_swap_nodes(ndord, nid, nextnid);
		}
		itnum ++;
	}while( isupdated & (itnum <itmax) );
	nds = ndord;
} 

