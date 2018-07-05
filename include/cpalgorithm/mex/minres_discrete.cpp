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
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif


//void bealgorithm( vector<int> indList, int N ){
std::vector<bool> minres_discrete( std::vector<int> indList, int N, int M ){
	//% --------------------------------
	//% Initialise
	//% --------------------------------
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
	
	std::vector<int> deg; deg.assign(N,0);
	for( auto ind : indList ){
		deg[ ind / N ] ++;
		deg[ ind % N ] ++;
	} 
	
	double Q = 0;
	//% init x
	if(N<=2){
		std::vector<bool> C; C.assign(N,true);
		return C;
	}
	std::vector<bool> C; C.assign(N,false);
	std::bernoulli_distribution rnd(0.5);
	int Ncore=0;
	for(int i = 0; i < N;i ++){
		int rand = rnd(mtrnd);
		if(rand==0){
			C[i] = false;
		}else{
			C[i] = true;
			Ncore ++;
		}
	} 
	
	//m = x'*A*x;
	std::vector<bool>fixed; fixed.assign(N,false);
	std::vector<bool>x = C;
	std::vector<bool>xbest; xbest.assign(N,false);
	int Ncorebest = 0; 
	//for j = 1:n	
	for( int j = 0;j < N;j++){
		//fixed =false(n,1);
		std::fill(fixed.begin(),fixed.end(),false);
		//Qold = -self.eval(G,x);
		//dQ = 0;dQmax = -Inf;
		double Mcore = 0;
		double Mperi = 0;
		for( auto ind : indList ){
			int cid = ind / N ;	
			int rid = ind % N ;
			if( C[rid] & C[cid] ) Mcore++;
			if( ~C[rid] & ~C[cid] ) Mperi++;
		}
		double Qold = -(Ncore * ( Ncore - 1)/2 - Mcore + Mperi); // pending
			
		double dQ = 0;
		double dQmax = -1 * std::numeric_limits<double>::max();
		int nid = 0;double qmax = 0;
		double q = 0;
		
		//for i = 1:ceil(n/2)
		for( int i = 0;i < ceil(N/2);i++){
			qmax = -1 * std::numeric_limits<double>::max();
			for(int k =0;k<N;k++){
				if( fixed[k] ) continue;
				q = -( (2*!!(x[k])-1) * (deg[k] -Ncore) + !!(x[k]) );
				if(qmax < q){
					nid = k;qmax = q;
				}
			}
			
			fixed[nid] = true;
			
			Ncore+= -(2*!!(x[nid])-1);
			
			x[nid] = !x[nid];
			
			dQ = dQ + qmax - Qold;
			Qold = qmax;
			
			if(dQmax < dQ){
				xbest = x;
				dQmax = dQ;
				Ncorebest = Ncore;
			}
		}
		if(dQmax <=0.00000001){
			break;
		}
		Ncore = Ncorebest;
		x = xbest;
		C = xbest;
	}
	return C;
} 



void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
   	int M = (int) mxGetPr(prhs[2])[0]; //length of edge list
	
	std::vector<int> L;
	for(int i = 0;i<M;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+M]-1);
		int ind = MIN(rid,cid) + MAX(rid,cid) * N;
		L.push_back(ind);
	}	
	std::vector<bool> C = minres_discrete(L,N,M);
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	double* R = mxGetPr(plhs[0]);
	for(int i = 0;i<N;i++){
		if( C[i]){
			R[i] = 1;
		}else{
			R[i] = 0;
		}	
	}	
}
