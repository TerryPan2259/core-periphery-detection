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

#include <math.h> 

class BEAlgorithm: public CPAlgorithm{
public:
	// Constructor 
	BEAlgorithm();
	BEAlgorithm(int num_runs);
	
	void detect(const Graph& G);
	
	void calc_Q(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<double>& x,
	    double& Q,
	    vector<double>& q);
	
protected: // function needed to be implemented
	int _num_runs; 
	void _detect_(const Graph& G, vector<double>& x, mt19937_64& mtrnd);
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
void BEAlgorithm::calc_Q(
    const Graph& G,
    const vector<int>& c,
    const vector<double>& x,
    double& Q,
    vector<double>& q)
{
	int N = G.get_num_nodes();
	double M = 0.0;
	double pa = 0;
	double pb = 0;
	double nc = 0;
	double mcc = 0;
	for( int i = 0;i < N;i++ ) {
		nc+=x[i];
	
		int sz = G.degree(i);	
		for( int k = 0; k < sz; k++ ) {
			Neighbour nei = G.get_kth_neighbour(i, k);
			int j = nei.get_node(); 
			double w = nei.get_w(); 
			mcc+=w * (x[i]+x[j] - x[i] * x[j]);
			M++;
		}
	}
	mcc = mcc/2;
	M = M /2;
	double M_b = (double)(nc * (nc-1) + 2 * nc * (N -nc))/2;
	pa = M / (double)(N * (N-1)/2);
	pb = M_b / (double)(N * (N-1)/2);	
	
	
	Q = ((double)mcc - pa * M_b ) / (sqrt(pa * (1-pa)) * sqrt(pb * (1-pb)) + 1e-30);
	Q = Q / (double)(N * (N-1)/2);
	
	if(Q > 1) Q = 1;
	if(Q < -1) Q = -1;

	vector<double> qtmp(1,Q);
	q= qtmp;
}

void BEAlgorithm::detect(const Graph& G){
    
    double Q = -1;
    int N = G.get_num_nodes();
    for (int i = 0; i < _num_runs; i++) {
        vector<int> ci(N, 0);
        vector<double> xi;
        vector<double> qi;
        double Qi = 0.0;
        _detect_(G, xi, _mtrnd);
		
        calc_Q(G, ci, xi, Qi, qi);
	
        if (Qi > Q) {
            _c = ci;
            _x = xi;
            _Q = Qi;
            _q = qi;
        }
    }
	
}

            
void BEAlgorithm::_detect_(const Graph& G, vector<double>& x, mt19937_64& mtrnd){
		
	// --------------------------------
	// Initialise _x randomly 
	// --------------------------------
	int N = G.get_num_nodes();
	double M = G.get_num_edges();
	double p = M / (double)(N * (N - 1) / 2 ); 
	
	vector<double> tmp(N, 0.0);
	x = tmp;
	uniform_real_distribution<double> dis(0.0, 1.0);	
	
	double Nperi = N;
	for(int i = 0;i<N; i++){
		if(dis(mtrnd) < 0.5) {
			x[i] = 1;
			Nperi-=1;
		}
	}
		
	// --------------------------------
	// Maximise the Borgatti-Everett quality function 
	// --------------------------------
	std::vector<double>xt = x;
	std::vector<double>xbest(N, 0.0);
	std::vector<bool>fixed(N, false);
	vector<double> Dperi(N, 0);

	for( int j = 0;j < N;j++){
		std::fill(fixed.begin(),fixed.end(),false);
		Nperi = 0.0;
		double numer = 0.0;
		for( int i = 0; i < N;i ++ ){
			Nperi+=(1-x[i]);
			Dperi[i] = 0;
			int sz = G.degree(i);
			for( int k = 0; k < sz;k ++ ){
				int nei = G.get_kth_neighbour(i, k).get_node();
				Dperi[i]+=1-x[nei];
				numer+= x[i] + x[nei] - x[i] * x[nei];
			}
		}

		numer = numer/2.0 -p*( (double)(N*(N-1.0))/2.0 - (double)Nperi*((double) Nperi-1.0)/2.0 );
		double pb = 1 -  (double)Nperi*(Nperi-1)/(double)(N*(N-1));
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
				double dnumer = (Dperi[k]- p * (Nperi-!!(1-xt[k])) ) * (2*(1-xt[k])-1);
				double newNperi = Nperi + 2*xt[k]-1;
				double pb = 1.0- (newNperi*(newNperi-1.0)) / (N*(N-1.0));
				double q = (numer + dnumer) / sqrt(pb*(1-pb));
				if( (qmax < q) & (pb*(1-pb)>0)){
					nid = k;qmax = q;numertmp = numer + dnumer;
				}
			}
			numer = numertmp;	
			Nperi+=2*xt[nid]-1;
	
			int sz = G.degree(nid);
			for(int k = 0;k<sz ;k++){
				int neik = G.get_kth_neighbour(nid, k).get_node();
				Dperi[ neik ]+=2*xt[nid]-1;
			}
		
			xt[ nid ] = 1-xt[ nid ];
			
			dQ = dQ + qmax - Qold;
			Qold = qmax;
	
			//% Save the core-periphery pair if it attains the largest quality	
			if(dQmax < dQ){
				xbest = xt;
				dQmax = dQ;
			}
			//fixed( nid ) = true; % Fix the label of node nid
			fixed[ nid ] = true; //% Fix the label of node nid
		}
		 	
		if (dQmax <= std::numeric_limits<double>::epsilon()){
			break;
		}
		
		xt = xbest; x = xbest;
	}
	x = xbest;
} 

