/*
*
* Header file of the BE algorithm
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

class BEAlgorithm: public BEAlgorithm{
public:
	// Constructor 
	BEAlgorithm();
	BEAlgorithm(int num_runs);
	
	void detect(const Graph& G);
	
	void calc_Q(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<bool>& x,
	    double& Q,
	    vector<double>& q);
	
protected: // function needed to be implemented
	int _num_runs; 
};



/*-----------------------------
Constructor
-----------------------------*/
BEAlgorithm::BEAlgorithm(int num_runs):CPAlgorithm(){
	BEAlgorithm();
	_num_runs = num_runs;
};

BEAlgorithm::BEAlgorithm(): CPAlgorithm(){
	_num_runs = 10;
};


/*-----------------------------
Functions inherited from the super class (CPAlgorithm)
-----------------------------*/
void BEAlgorithm::detect(const Graph& G){
	_km_modmat_louvain(G.to_matrix(), _num_runs, _c, _x, _Q, _q, _mtrnd);
}

void BEAlgorithm::calc_Q(
    const Graph& G,
    const vector<int>& c,
    const vector<bool>& x,
    double& Q,
    vector<double>& q)
{
	vector<vector<double>>M = G.to_matrix();
	_calc_Q_modmat(M,c,x,Q,q);
}

/*-----------------------------
Private functions (internal use only)
-----------------------------*/
