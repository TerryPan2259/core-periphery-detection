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
	void _detect_(const Graph& G, vector<double>& x, mt19937_64& mtrnd);
	void _community_detection(const Graph& G, vector<int>& c, vector<double>& x, mt19937_64& mtrnd);
	void _louvain(const Graph& G, const int num_of_runs, vector<vector<bool>>& xlist, mt19937_64& mtrnd);
	void _louvain_core(const Graph& G, vector<int>& C, const double M, mt19937_64& mtrnd);
	double _calc_dQmod( double d, double degi, double D, double selfw, const double M );
	double _calc_Qmod(const Graph& G, vector<int>& C, const double M);
	void _coarsing(const Graph& G, const vector<int>& c,Graph& newG);
	void _modularity_label_switching(const Graph& G,vector<int>& C, const double M,mt19937_64& mtrnd);
	
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
	int N = G.get_num_nodes();
	double M = 0.0;
	double pa = 0;
	double pb = 0;
	int nc = 0;
	int mcc = 0;
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
	
	
	Q = ((double)mcc - pa * M_b ) / (sqrt(pa * (1-pa)) * sqrt(pa * (1-pb)));
	Q = Q / (double)(N * (N-1)/2);
	vector<double> qtmp(1,Q);
	q= qtmp;
}

void Divisive::detect(const Graph& G){
    
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

            
void Divisive::_detect_(const Graph& G, vector<double>& x, mt19937_64& mtrnd){
		
	// --------------------------------
	// Initialise _x randomly 
	// --------------------------------
	int N = G.get_num_nodes();
	double M = G.get_num_edges();
	double p = M / (double)(N * (N - 1) / 2 ); 
	
	vector<double> tmp(N, 0.0);
	x = tmp;
	uniform_real_distribution<double> dis(0.0, 1.0);	
	
	int Nperi = N;
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
	vector<int> Dperi(N, 0);

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


double Divisive::_calc_dQmod( double d, double degi, double D, double selfw, const double M ){
	return 2*( d - degi*D/(2.0*M)) + (selfw-degi*degi/(2.0*M)); 
}

double Divisive::_calc_Qmod(
	const Graph& G, 
	vector<int>& C, 
	const double M
	 ){
	
	double retval = 0;
	int N = C.size();
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
		
        int N = c.size();
    	int maxid = 0;
	
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
	
        int N=C.size();
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
			
			double dQold = calc_dQmod( toC[ncid], deg[nid], D[ncid] - deg[nid], selfw, M );
			double dQ = 0;
			for(int j = 0;j<neighbourNum;j++){
				Neighbour nb = G.get_kth_neighbour(nid, j);
				int nei = nb.get_node();
				double w = nb.get_w();
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
 		int labsize = labs.size();
		for(int j=0;j<labsize;j++){
			if(labs[j]==C[i]){
				cid = j;
				break;
			}
		}
		
		if (cid<0) {
			labs.push_back(C[i]);
			cid = labs.size()-1;
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
	unsigned int prevGraphSize = C.size();
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
		int Ctsize = Ct.size();
		for(int i = 0;i<Ctsize;i++){
			Ct[i] = Zt[ Ct[i] ];
		}
		
		if(Qt>Qbest){
			C = Ct;
			Qbest = Qt;
		}
	}while( newG.get_num_nodes()!= prevGraphSize);
		
}

void Divisive::_louvain(
    const Graph& G,
    const int num_of_runs,
    vector<vector<bool>>& xlist,
    mt19937_64& mtrnd
    ){

    int N = G.get_num_nodes();
    vector<int> C(N);
    int K = xlist.size();
    for(int k = 0;k < K;k++){
        xlist[k].clear();
    }
    xlist.clear();

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
    K = *max_element(cbest.begin(),cbest.end()) + 1; 
    for(int k = 0; k < K; k++){
	vector<bool> tmp(N);
    	for(int i = 0; i < N; i++){
		tmp[i] = cbest[i]==k;
    	}	
	xlist.push_back(tmp);	
    }
}
