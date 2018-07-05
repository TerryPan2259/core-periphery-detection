#include <array>
#include <random>
#include <algorithm>
#include <bitset>
#include <vector>
#include <set>
#include "mex.h"
#include "math.h"
#include <random>
#include <iostream>     // std::cout, std::end
#include <fstream>
#include <limits>
#include <cmath>
#include <cfloat>
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

void bealgorithm(vector<std::vector<int>>& adjList,vector<bool>& X, const int N, const double M, const double p){
	// --------------------------------
	// Initialise
	// --------------------------------
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
	
		
	// --------------------------------
	// Initialise X using MINRES algorithm 
	// --------------------------------
	std::vector<int> deg; deg.assign(N,0);
	for( int i = 0;i < N;i++ ) deg[i] = adjList[i].size();
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



void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
   	double M = (double) mxGetPr(prhs[2])[0]; //length of edge list
	vector<vector<int>> adjList;	
	for(int i = 0;i < N;i++){
		std::vector<int> tmp;
		adjList.push_back(tmp);
	}
	for(int i = 0;i<M;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+(int)M]-1);
		if(rid==cid){
		}else{
			adjList[rid].push_back(cid);
			adjList[cid].push_back(rid);
		}
	}
	double p = M / ((double)N*((double)N-1.0)/2.0);
	vector<bool> X;X.assign(N,false);
	bealgorithm(adjList, X, N, M, p );
    	//plhs[0] = mxXreateLogicalMatrix( (mwSize)N, (mwSize)1); 
    	//bool* R = mxGetLogicals(plhs[0]);
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	double* R = mxGetPr(plhs[0]);
	for(int i = 0;i<N;i++) R[i] = !!(X[i]);
	
	for(int i = 0;i < N;i++){
		adjList[i].clear();
	}
}
