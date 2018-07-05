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

std::vector<std::vector<int>> adjList;
std::vector<std::vector<double>> wList;
std::vector<bool> C;
std::vector<bool> Cbest;
std::vector<int> P;
int N;
int M;
double p;
double beta = 10;
double calcQ(){
	double retval = 0;
	for(int i = 0;i < N;i++){
		if(C[i]==true) continue;
		for(int j = 0;j<adjList[i].size();j++){
			int nei = adjList[i][j];
			if(nei!=P[i]) continue;
			retval+=wList[i][j] - p;
		}
	}
	return retval;	
}

void findPeriphery(){
	for(int i = 0;i < N;i++){
		if(C[i]==true) {
			P[i] = i;
			continue;
		}
		int core = -1;double w = -1;
		for(int j = 0;j<adjList[i].size();j++){
			int nei = adjList[i][j];
			if( C[nei] & w < wList[i][j] ){
				w = wList[i][j];	
				core = nei;
			}
		}
		if (w-p >0.0){
			P[i] = core;
		}else{
			P[i] = -1;
		}
	}
} 

//void bealgorithm( vector<int> indList, int N ){
void starkm_kl(int maxStableLoopNum,int maxItNum){
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
	std::uniform_real_distribution<double> udist(0.0,1.0);
	std::uniform_int_distribution<int> idist(0,N-1);
	
	// ===========================
	// Detect core-periphery pairs using a MCMC algorithm 
	// ===========================
	// initialise core-periphery pairs
	C.assign(N,false);P.assign(N,0);	
	for(int i = 0;i < N;i++){
		C[i] = udist(mtrnd)<=0.5;
	}
	findPeriphery();
	double Qmax = calcQ();
	// initialise parameters
	
	std::vector<bool> fixed;
	std::vector<bool> Ctmp;
	for(int i = 0;i<N;i++){
		fixed.assign(N,false);
		C = Cbest;
		findPeriphery();
		double Qold = calcQ();
		double qmax = Qold;
		Ctmp = C;
		for(int j = 0;j<N;j++){
			double q = 0;
			int maxid = -1;
			for(int k = 0;k<N;k++){
				if(fixed[k]) continue;
				C[k] = !C[k];
				findPeriphery();
				double qk = calcQ();
				if(q < qk){
					q = qk;
					maxid = k;
				}
				C[k] = !C[k];
			}
			
			if(maxid<0) continue;	
			fixed[maxid] = true;	
			C[maxid] = !C[maxid];
			
			if(qmax < q){
				Ctmp = C;
				qmax = q;
			}
		}
		if(qmax <= Qmax) break;
		Cbest = Ctmp;
		Qmax = qmax;
	}	
	
	C = Cbest;
	findPeriphery();
} 



void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	N = (int) mxGetPr(prhs[1])[0];
   	M = (int) mxGetPr(prhs[2])[0]; // length of edge list
   	int maxItNum = (int) mxGetPr(prhs[3])[0]; // maximum number of iterations
   	int maxStableLoopNum = (int) mxGetPr(prhs[4])[0]; // maximum number of intervals for updating Qmax 
	for(int i = 0;i < N;i++){
		std::vector<int> tmp;
		adjList.push_back(tmp);
		std::vector<double> tmp2;
		wList.push_back(tmp2);
	}
	p = 0;
	for(int i = 0;i<M;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+M]-1);
		double w = edgeList[i+2*M];
		adjList[rid].push_back(cid);
		adjList[cid].push_back(rid);
		wList[rid].push_back(w);
		wList[cid].push_back(w);
		p = p + w;
	}
	p = 2.0 * p / (double)((N*(N-1)));
		
	starkm_kl(maxStableLoopNum,maxItNum);
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	plhs[1] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
	
    	double* retC = mxGetPr(plhs[0]);
    	double* retX = mxGetPr(plhs[1]);
	
	std::vector<int> labs;
	for(int i=0;i<N;i++){
		int cid = -1;
		for(int j=0;j<labs.size();j++){
			if(labs[j]==P[i]){
				cid = j+1;
				break;
			}
		}
		
		if (cid<0) {
			labs.push_back(P[i]);
			cid = labs.size();
		}
		P[i] = cid;		
	}
	for(int i = 0;i<N;i++){
		retC[i] = P[i];
		retX[i] = !!(C[i]);	
	}
	for(int i = 0;i < N;i++){
		adjList[i].clear();
		wList[i].clear();
	}
	adjList.clear();	
	wList.clear();	
	C.clear();	
	Cbest.clear();	
	P.clear();	
}
