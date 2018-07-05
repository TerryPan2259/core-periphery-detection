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

class MINRES: public MINRES{
public:
	// Constructor 
	MINRES();
	
	void detect(const Graph& G);
	
	void calc_Q(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<bool>& x,
	    double& Q,
	    vector<double>& q);
private:

	vector<int> MINRES::_sortIndex(const vector<int>& Qs);
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
	vector<bool> tmp2(N, false);
	for(int k = 0;k<=kbest;k++){
		tmp2[ ord[k] ] = true;
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
    const vector<bool>& x,
    double& Q,
    vector<double>& q)
{
	
	Q = 0.0;
	int mcc=0;
	int mpp = 0;
	int ncc = 0;

	int N = G.get_num_nodes();
	for(int i = 0; i < N; i++){
		int sz = G.degree(i);
		for(int j = 0; j < sz; j++){
			int nei = -1; double w = -1;
			G.get_weight(i, j, nei, w);
			if(x[i] & x[j]){	
				mcc++;
			}
			if(!x[i] & !x[j]){	
				mpp++;
			}
		}
		if(x[i]) ncc+=1;
	}
	
	//(ncc * ncc - mcc) number of absent edges in core
	//(mpp) number of present edges in periphery 
	Q = (ncc * ncc - mcc) + mpp;
	Q = -Q;
	q = Q;	
}

/*-----------------------------
Private functions (internal use only)
-----------------------------*/
vector<int> MINRES::sortIndex(const vector<int>& Qs){
    vector<int> y(Qs.size());
    size_t n(0);
    generate(std::begin(y), std::end(y), [&]{ return n++; });

    sort(  std::begin(y), 
                std::end(y),
                [&](int i1, int i2) { return Qs[i1] > Qs[i2]; } );
    return y;
}


