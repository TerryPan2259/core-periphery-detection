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
