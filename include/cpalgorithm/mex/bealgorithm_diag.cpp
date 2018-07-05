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
std::vector<bool> bealgorithm_diag( std::vector<int> indList, int N, int M ){
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
	M = indList.size();
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
	std::vector<int> Dcore; Dcore.assign( N, 0 );
	int Mcc = 0; int Mpp = 0;
	for( int ind : indList ){
		int cid = ind / N ;	
		int rid = ind % N ;
		if( C[cid] & C[rid] ) Mcc++;
		if( !C[cid] & !C[rid] ) Mpp++;
		
		if( C[cid] ){
			Dcore[rid]++;
		}	
		if( C[rid] ){
			Dcore[cid]++;
		}	
	}
	std::vector<bool>x = C;
	std::vector<bool>xbest; xbest.assign(N,false);
	std::vector<bool>fixed; fixed.assign(N,false);
	std::vector<int>Dcorebest; Dcorebest.assign(N,0);
	int Ncorebest = 0;
	int Mccbest = Mcc; int Mppbest = Mpp;
	// --------------------------------
	// Maximise the Borgatti-Everett quality function 
	// --------------------------------
	for( int j = 0;j < N;j++){
		std::fill(fixed.begin(),fixed.end(),false);
		int nid = 0;double qmax, q, dQ=0;
		double dQmax = -1 * std::numeric_limits<double>::max();
		double pa = (Mcc + Mpp)/(Ncore*(Ncore-1)/2  + ( N - Ncore )*( N - Ncore -1)/2);
		double pb = (Ncore*(Ncore-1)/2)/(Ncore*(Ncore-1)/2  + ( N - Ncore )*(N - Ncore -1)/2);
		double Qold = (1-pa)*Mcc - ((Ncore-1)*Ncore/2-Mcc)*pa;
		Qold = Qold /( sqrt(pa*(1-pa)*pb*(1-pb)) );
		if(isinf(Qold)) Qold = 0;
		
		for( int i = 0;i < ceil(N/2);i++){
			double qmax = -1 * std::numeric_limits<double>::max();
			for(int k = 0; k < N;k++){
				if(fixed[k]) continue;
				double nc = Ncore +(1-2*!!(x[k]));
				double np = N-nc;
				int Ec = Mcc + (1-2*!!(x[k]))*Dcore[k];
				int Ep = Mpp + (1-2*( 1 - !!(x[k]) ) )*(deg[k]-Dcore[k]);
				pa = ( Ep + Ec ) / ( nc * ( nc-1 ) / 2  + np * ( np - 1 ) / 2 );
				pb = (nc*(nc-1)/2)/(nc*(nc-1)/2  + np*(np-1)/2);
				double q = (1-pa)*Ec - ((nc-1)*nc/2-Ec)*pa;
				q = q / ( sqrt(pa*(1-pa)*pb*(1-pb)) );
				if(qmax < q & ~isinf( q )){
					nid = k;qmax = q;
				}
			}
				
			Mcc = Mcc + ( 1 - 2 *!!(x[nid]) )*Dcore[nid];
			Mpp = Mpp + (1-2*( 1 - !!(x[nid]) ))*(deg[nid]-Dcore[nid]);
			for(auto ind : indList ){
				int rid = ind / N;
				int cid = ind % N;
				if(rid==nid) Dcore[cid]+=-2 * !!(x[nid]) + 1;
				if(cid==nid) Dcore[rid]+=-2 * !!(x[nid]) + 1;
			} 
			Ncore = Ncore + (1-2*!!(x[nid]));
			x[ nid ] = !x[ nid ];
			
			dQ = dQ + qmax - Qold;
			Qold = qmax;
			//% Save the core-periphery pair if it attains the largest quality	
			if(dQmax < dQ){
				xbest = x;
				dQmax = dQ;
				Ncorebest = Ncore; Dcorebest = Dcore;
				Mccbest = Mcc;Mppbest = Mpp;
			}
			fixed[ nid ] = true; //% Fix the label of node nid
		}
		 	
		if (dQmax <= std::numeric_limits<double>::epsilon()){
			break;
		}
		Ncore = Ncorebest; Dcore = Dcorebest;
		Mcc = Mccbest;Mpp = Mppbest;
		x = xbest; C = xbest;
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
	std::vector<bool> C = bealgorithm_diag(L,N,M);
    	//plhs[0] = mxCreateLogicalMatrix( (mwSize)N, (mwSize)1); 
    	//bool* R = mxGetLogicals(plhs[0]);
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	double* R = mxGetPr(plhs[0]);
	for(int i = 0;i<N;i++){
		if( C[i]){
			R[i] = 1;
		}else{
			R[i] = 0;
		}	
		//R[i] = C[i];
	}	
}
