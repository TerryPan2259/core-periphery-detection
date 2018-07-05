#include <array>
#include <random>
#include <algorithm>
#include <vector>
#include "math.h"
#include <random>
#include <iostream>     // std::cout, std::end
#include <fstream>
#include <limits>
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif
using namespace std;

vector<int> sortIndex(const vector<int>& Qs){
    vector<int> y(Qs.size());
    size_t n(0);
    generate(std::begin(y), std::end(y), [&]{ return n++; });

    sort(  std::begin(y), 
                std::end(y),
                [&](int i1, int i2) { return Qs[i1] > Qs[i2]; } );
    return y;
}

void degree_algorithm(vector<std::vector<int>>& adjList,vector<bool>& X, const vector<int>& deg,const int N, const double M){
	double mu = 2*M/(double)N;	
	for(int k = 0;k<=N;k++){
		if( mu < deg[k]){
			X[ k ] = true;
		}else{
			X[ k ] = false;
		}
	}
}

void lip_algorithm(vector<std::vector<int>>& adjList,vector<bool>& X, const vector<int>& deg,const int N, const double M){
	vector<int> ord = sortIndex(deg);	
	double Z = M;double Zbest = numeric_limits<double>::max();
	int kbest = 0;
	for(int k = 0;k<N;k++){
		Z = Z + k - 1 - deg[ ord[k] ];
		if(Z < Zbest){
			kbest = k;
			Zbest = Z;
		}
		X[k]=false;
	}
	
	for(int k = 0;k<=kbest;k++){
		X[ ord[k] ] = true;
	}
	int Nperi = N - kbest -1;
}

double calcQbe(vector<vector<int>>& adjList, vector<bool>& X, const int N, const double M, const double p){
	double numer = 0;
	int ncore = 0;
	for(int i = 0;i < adjList.size();i++){
		if(X[i]) ncore++;
	}
	double pcp = (double)(N*(N-1) - (N-ncore) * ((N-ncore)-1))/(double)(N*(N-1));
	for(int i = 0;i < adjList.size();i++){
		for(int j = 0;j < adjList[i].size();j++){
			if(X[i] | X[ adjList[i][j] ] ) {
				numer+=1-pcp;	
			}else{
				numer+=-pcp;
			}
		}
	}
	double denom = sqrt( (double)N*((double)N-1)*p*(1-p) ) * sqrt( (double)N*((double)N-1)*pcp*(1-pcp)  );
	double Q = 0;
	if(denom<1e-20){
		Q = 0;
	}else{
		Q = numer / denom;
	}
	return Q;	
}

double calcQbe_ignore_cp_edges(vector<vector<int>>& adjList, vector<bool>& X, const int N, const double M, const double p){
	double numer = 0;
	int ncore = 0;
	for(int i = 0;i < adjList.size();i++){
		if(X[i]) ncore++;
	}
	double pcp = (double)(ncore * (ncore-1))/(double)(N*(N-1));
	for(int i = 0;i < adjList.size();i++){
		for(int j = 0;j < adjList[i].size();j++){
			if(X[i] & X[ adjList[i][j] ] ) {
				numer+=1-pcp;	
			}else{
				numer+=-pcp;
			}
		}
	}
	double denom = sqrt( (double)N*((double)N-1)*p*(1-p) ) * sqrt( (double)N*((double)N-1)*pcp*(1-pcp)  );
	double Q = 0;
	if(denom<1e-20){
		Q = 0;
	}else{
		Q = numer / denom;
	}
	return Q;	
}

void gen_rand_er(vector<vector<int>>& adjList,const int N, const int M, std::mt19937_64& mtrnd) {
			
	for( int mid = 0;mid<M;mid++){
		int rid=-1;
		int cid=-1;
		while(true){	
			rid = std::uniform_int_distribution<>(0, N-1)(mtrnd);
			cid = std::uniform_int_distribution<>(0, N-1)(mtrnd);
			if(rid >cid){
				int tmp = rid;
				rid = cid;
				cid = tmp;
			}
		
			if(rid == cid) continue;
			bool skip = false;
			for(int i =0;i < adjList[rid].size();i++){
				if(adjList[rid][i]==cid){
					skip = true;
					break;
				}
			}
			if( skip ) continue;
			break;
		}
		
		adjList[rid].push_back(cid);
		adjList[cid].push_back(rid);
		
	};
}

void bealgorithm(vector<std::vector<int>>& adjList,vector<bool>& X, const int N, const double M, const double p, std::mt19937_64& mtrnd ){	
		
	// --------------------------------
	// Initialise X using MINRES algorithm 
	// --------------------------------
	std::vector<int> deg; deg.assign(N,0);
	for( int i = 0;i < N;i++ ) deg[i] = adjList[i].size();
	lip_algorithm(adjList,X,deg,N,M);
	double Nperi = 0;
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
			for( int k = 0; k < adjList[i].size();k ++ ){
				int nei = adjList[i][k];
				if( !X[nei] ) Dperi[i]++;
				if(X[i] | X[adjList[i][k]]) numer++;	
			}
		}
		numer = numer/2.0 -p*( N*(N-1.0)/2.0 - (double)Nperi*((double) Nperi-1.0)/2.0 );
		double pb = 1 -  Nperi*(Nperi-1)/(N*(N-1));
		double Qold = numer / sqrt(pb*(1-pb));
		
		double dQ = 0;
		double dQmax = -1 * std::numeric_limits<double>::max();
		int nid = 0;
		
		for( int i = 0;i < N/2;i++){
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
			for(int k = 0;k<adjList[nid].size();k++){
				Dperi[ adjList[nid][k] ]+=2*!!(x[nid])-1;
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

void bealgorithm_without_cp_edges(vector<std::vector<int>>& adjList,vector<bool>& X, const int N, const double M, const double p, std::mt19937_64& mtrnd ){
	
	// --------------------------------
	// Initialise X using MINRES algorithm 
	// --------------------------------
	std::vector<int> deg; deg.assign(N,0);
	for( int i = 0;i < N;i++ ) deg[i] = adjList[i].size();
	lip_algorithm(adjList,X,deg,N,M);
	double Ncore = 0;
	// --------------------------------
	// Maximise the Borgatti-Everett quality function 
	// --------------------------------
	std::vector<bool>x = X;
	std::vector<bool>xbest; xbest.assign(N,false);
	std::vector<bool>fixed; fixed.assign(N,false);
	vector<int> toC; toC.assign(N,0);
	for( int j = 0;j < N;j++){
		std::fill(fixed.begin(),fixed.end(),false);
		Ncore = 0.0;
		double Qcc = 0.0; // number of absent edges within periphery
		for( int i = 0; i < N;i ++ ){
			if(X[i]) Ncore++;
			toC[i] = 0;
			for( int k = 0; k < adjList[i].size();k ++ ){
				int nei = adjList[i][k];
				if(X[nei]) toC[i]++;
				if(X[i] & X[adjList[i][k]]) Qcc++;	
			}
		}
		Qcc = Qcc/2.0 - p*( (double)Ncore*((double) Ncore-1.0)/2.0 );
		double pb = Ncore*(Ncore-1)/(N*(N-1));
		double Qold = Qcc / sqrt(pb*(1-pb));
		
		double dQ = 0;
		double dQmax = -1 * std::numeric_limits<double>::max();
		int nid = 0;
		
		for( int i = 0;i < N/2;i++){
			double qmax = -1 * std::numeric_limits<double>::max();
			
			// select a node of which we update the label 
			double Qcctmp = Qcc;
			for(int k =0;k<N;k++){
				if( fixed[k] ) continue;
				double dQcc = (toC[k]- p * (Ncore-!!(x[k])) ) * (2*!!(!x[k])-1); // number of edges incremented after swithing the class of node k
			
				double newNcore = Ncore + 2*(!!(!x[k]))-1;
					
				double pb = (newNcore*(newNcore-1.0)) / (N*(N-1.0));
				double q = (Qcc + dQcc) / sqrt(pb*(1-pb));
				if(qmax < q & pb*(1-pb)>0){
					nid = k;qmax = q;Qcctmp = Qcc + dQcc;
				}
			}
			
			Qcc = Qcctmp;	
			Ncore+=2*!!(!x[nid])-1;
			for(int k = 0;k<adjList[nid].size();k++){
				toC[ adjList[nid][k] ]+=2*!!(!x[nid])-1;
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

double estpval_boyd_stest(vector<vector<int>>& adjList,vector<bool>& X, const int N, const double M, const double p,const int numRandGraph,double pval){
	/*
	 * Initialise	
	 * 
	 */
	std::mt19937_64 mtrnd;
    	int seeds[624];
       	size_t size = 624*4; //Declare size of data
       	std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary); //Open stream
       	if (urandom) //Check if stream is open
       	{
       	    urandom.read(reinterpret_cast<char*>(seeds), size); //Read from urandom
       	    urandom.close(); //close stream
       	}
       	else //Open failed
       	{
            		std::cerr << "Failed to open /dev/urandom" << std::endl;
       	}
    	std::seed_seq seed(&seeds[0], &seeds[624]);
    	mtrnd.seed(seed);

	vector<vector<int>> randGraph;	
	
	double Qth = calcQbe(adjList, X, N, M, p);
	vector<bool> Xrand(N);Xrand.assign(N,false);
	double estpval = 0;
	for(int i = 0;i < numRandGraph;i++){
		vector<vector<int>>().swap(randGraph);
		for(int i = 0;i < N;i++){
			std::vector<int> tmp;
			randGraph.push_back(tmp);
		}
		gen_rand_er(randGraph,N, M, mtrnd);
		bealgorithm(randGraph,Xrand,N,M,p,mtrnd);
		double Qrand = calcQbe(randGraph,Xrand,N,M,p);
		if(Qth <=Qrand){
			estpval+= (double)1.0/(double)numRandGraph; 
		}
		if(estpval > pval){
			estpval = estpval * numRandGraph / (i+1); 
			return estpval;
		}
	}
	return estpval;
}


double estpval_boyd_stest_for_dense_net(vector<vector<int>>& adjList,vector<bool>& X, const int N, const double M, const double p,const int numRandGraph,double pval){
	/*
	 * Initialise	
	 * 
	 */
	std::mt19937_64 mtrnd;
    	int seeds[624];
       	size_t size = 624*4; //Declare size of data
       	std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary); //Open stream
       	if (urandom) //Check if stream is open
       	{
       	    urandom.read(reinterpret_cast<char*>(seeds), size); //Read from urandom
       	    urandom.close(); //close stream
       	}
       	else //Open failed
       	{
            		std::cerr << "Failed to open /dev/urandom" << std::endl;
       	}
    	std::seed_seq seed(&seeds[0], &seeds[624]);
    	mtrnd.seed(seed);

	vector<vector<int>> randGraph;	
	double Qth = calcQbe(adjList, X, N, M, p);
	vector<bool> Xrand(N);Xrand.assign(N,false);
	double estpval = 0;
	for(int i = 0;i < numRandGraph;i++){
		
		vector<vector<int>>().swap(randGraph);
		for(int i = 0;i < N;i++){
			std::vector<int> tmp;
			randGraph.push_back(tmp);
		}
		gen_rand_er(randGraph,N, N*(N-1)/2-M, mtrnd); // generate complementary graphs
		bealgorithm_without_cp_edges(randGraph,Xrand,N,N*(N-1)/2-M,1-p,mtrnd); // detect periphery subsets
		double Qrand = calcQbe_ignore_cp_edges(randGraph,Xrand,N,N*(N-1)/2-M,1-p);
		if(Qth <=Qrand){
			estpval+= (double)1.0/(double)numRandGraph; 
		}
		if(estpval > pval){
			estpval = estpval * numRandGraph / (i+1); 
			return estpval;
		}
	}
	return estpval;
}
