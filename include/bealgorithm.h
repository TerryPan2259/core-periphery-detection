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
	
	int N = G.get_num_nodes();
	int M = G.get_num_edges();
	
	double p = (double) M / (double)(N * (N-1) / 2)
	vector<int> tmp(N,0);	
	
	_detect_(G, _x, N, M, p, _mtrnd);
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
vector<int> sortIndex(const vector<int>& Qs){
    vector<int> y(Qs.size());
    size_t n(0);
    generate(std::begin(y), std::end(y), [&]{ return n++; });

    sort(  std::begin(y), 
                std::end(y),
                [&](int i1, int i2) { return Qs[i1] > Qs[i2]; } );
    return y;
}

void BEAlgorithm::_detect_(Graph& G,vector<bool>& X, const int N, const double M, const double p, mt19937_64 mtrnd){
//void bealgorithm(vector<std::vector<int>>& adjList,vector<bool>& X, const int N, const double M, const double p, mt19937_64 mtrnd){
		
	// --------------------------------
	// Initialise X using MINRES algorithm 
	// --------------------------------
	std::vector<int> deg; deg.assign(N,0);
	for( int i = 0;i < N;i++ ) deg[i] = G.degree(i);

	vector<int> ord = sortIndex(deg);	
	double Z = M;double Zbest = numeric_limits<double>::max();
	int kbest = 0;
	for(int k = 0;k<N;k++){
		Z = Z + k - 1 - deg[ ord[k] ];
		if(Z < Zbest){
			kbest = k;
			Zbest = Z;
		}
	}
	
	X.assign(N,false);
	for(int k = 0;k<=kbest;k++){
		X[ ord[k] ] = true;
	}
	int Nperi = N - kbest -1;
	
	// --------------------------------
	// Maximise the Borgatti-Everett quality function 
	// --------------------------------
	std::vector<bool>x = X;
	std::vector<bool>xbest; xbest.assign(N,false);
	std::vector<bool>fixed; fixed.assign(N,false);
	vector<int> Dperi; Dperi.assign(N,0);
	for( int j = 0;j < N;j++){
		std::fill(fixed.begin(),fixed.end(),false);
		Nperi = 0.0;
		double numer = 0.0;
		for( int i = 0; i < N;i ++ ){
			if(!X[i]) Nperi++;
			Dperi[i] = 0;
			int sz = G.degree(i);
			for( int k = 0; k < sz;k ++ ){
				double w;
				G.get_weight(i, k, nei, w)
				if( !X[nei] ) Dperi[i]++;
				if(X[i] | X[nei]) numer++;	
			}
		}
		numer = numer/2.0 -p*( N*(N-1.0)/2.0 - (double)Nperi*((double) Nperi-1.0)/2.0 );
		double pb = 1 -  Nperi*(Nperi-1)/(N*(N-1));
		double Qold = numer / sqrt(pb*(1-pb));
		
		
		double dQ = 0;
		double dQmax = -1 * std::numeric_limits<double>::max();
		int nid = 0;
		
		for( int i = 0;i < N;i++){
			double qmax = -1 * std::numeric_limits<double>::max();
			
			// select a node of which we update the label 
			double numertmp = numer;
			for(int k =0;k<N;k++){
				if( fixed[k] ) continue;
				double dnumer = (Dperi[k]- p * (Nperi-!!(!x[k])) ) * (2*!!(!x[k])-1);
				double newNperi = Nperi + 2*(!!x[k])-1;
				double pb = 1.0- (newNperi*(newNperi-1.0)) / (N*(N-1.0));
				double q = (numer + dnumer) / sqrt(pb*(1-pb));
				if(qmax < q & pb*(1-pb)>0){
					nid = k;qmax = q;numertmp = numer + dnumer;
				}
			}
			numer = numertmp;	
			Nperi+=2*!!(x[nid])-1;
	
			int sz = G.degree(i);
			for(int k = 0;k<sz ;k++){
				int neik = 0; double wk;
				G.get_weight(nid, k, neik, wk);
				Dperi[ neik ]+=2*!!(x[nid])-1;
			}
		
			x[ nid ] = !x[ nid ];
			
			dQ = dQ + qmax - Qold;
			Qold = qmax;
	
			//% Save the core-periphery pair if it attains the largest quality	
			if(dQmax < dQ){
				xbest = x;
				dQmax = dQ;
			}
			//fixed( nid ) = true; % Fix the label of node nid
			fixed[ nid ] = true; //% Fix the label of node nid
		}
		 	
		if (dQmax <= std::numeric_limits<double>::epsilon()){
			break;
		}
		
		x = xbest; X = xbest;
	}
	X = xbest;
} 

