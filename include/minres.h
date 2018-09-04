/*
*
* Header file of the MINRES algorithm
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

class MINRES: public CPAlgorithm{
public:
	// Constructor 
	MINRES();
	
	void detect(const Graph& G);
	
	void calc_Q(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<double>& x,
	    double& Q,
	    vector<double>& q);
private:
        vector<int> _sortIndex(const vector<int>& Qs);
};



/*-----------------------------
Constructor
-----------------------------*/
MINRES::MINRES():CPAlgorithm(){
};



/*-----------------------------
Functions inherited from the super class (CPAlgorithm)
-----------------------------*/
void MINRES::detect(const Graph& G){
	
	int N = G.get_num_nodes();
	double M = 0.0;
	std::vector<int> deg(N, 0);
	for( int i = 0;i < N;i++ ) {
		deg[i] = G.degree(i);
		M +=deg[i];
	}

	vector<int> ord = _sortIndex(deg);
	double Z = M;
	double Zbest = numeric_limits<double>::max();
	int kbest = 0;
	for(int k = 0;k<N;k++){
		Z = Z + k - 1 - deg[ ord[k] ];
		if(Z < Zbest){
			kbest = k;
			Zbest = Z;
		}
	}
	
	// set
	vector<double> tmp2(N, 0.0);
	for(int k = 0;k<=kbest;k++){
		tmp2[ ord[k] ] = 1.0;
	}
	_x = tmp2;
	
	// set
	vector<int> tmp(N,0);
	_c = tmp;	
	
	calc_Q(G, _c, _x, _Q, _q);
}

void MINRES::calc_Q(
    const Graph& G,
    const vector<int>& c,
    const vector<double>& x,
    double& Q,
    vector<double>& q)
{
	
	Q = 0.0;
	double mcc=0;
	double mpp = 0;
	double ncc = 0;

	int N = G.get_num_nodes();
	for(int i = 0; i < N; i++){
		int sz = G.degree(i);
		for(int j = 0; j < sz; j++){
			int nei = -1; double w = -1;
			G.get_weight(i, j, nei, w);
			mcc+=x[i] * x[j];
			mpp+= (1-x[i]) * (1- x[j]);
		}
		ncc+=x[i];
	}
	
	//(ncc * ncc - mcc) number of absent edges in core
	//(mpp) number of present edges in periphery 
	Q = (ncc * ncc - mcc) + mpp;
	Q = -Q;
	
	vector<double> tmp(1, Q);
	q = tmp;
}

/*-----------------------------
Private functions (internal use only)
-----------------------------*/
vector<int> MINRES::_sortIndex(const vector<int>& data){
	vector<int> index((int)data.size(), 0);
	for (int i = 0 ; i != index.size() ; i++) {
    		index[i] = i;
	}
	sort(index.begin(), index.end(),
    		[&](const int& a, const int& b) {
        	return (data[a] > data[b]);
    		}
	);
    return index;
}
