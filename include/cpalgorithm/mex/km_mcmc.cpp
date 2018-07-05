#include "mex.h"
#include "km_mcmc.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
   	int M = (int) mxGetPr(prhs[2])[0]; // length of edge list
   	int K = (int) mxGetPr(prhs[3])[0]; // maximum number of core-periphery pairs
   	int maxItNum = (int) mxGetPr(prhs[4])[0]; // maximum number of iterations
   	int maxStableLoopNum = (int) mxGetPr(prhs[5])[0]; // maximum number of intervals for updating Qmax 
	int maxNc = (int) mxGetPr(prhs[6])[0]; // maximum number of core nodes in a core-periphery pair
	double beta  = (double) mxGetPr(prhs[7])[0];
	
	
	vector<vector<int>> adjList;
	vector<vector<double>> wList;
	for(int i = 0;i < N;i++){
		std::vector<int> tmp;
		adjList.push_back(tmp);
		std::vector<double> tmp2;
		wList.push_back(tmp2);
	}
	double p = 0;
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
	p = 2 * p / (N*(N-1));
		
	
	
	vector<int> C;vector<bool> X;
	km_mcmc(adjList, wList, p, C, X, K, maxNc, beta, maxStableLoopNum,maxItNum);
	
	
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	plhs[1] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
	
    	double* retC = mxGetPr(plhs[0]);
    	double* retX = mxGetPr(plhs[1]);
	
	std::vector<int> labs;
	for(int i=0;i<N;i++){
		int cid = -1;
		for(int j=0;j<labs.size();j++){
			if(labs[j]==C[i]){
				cid = j+1;
				break;
			}
		}
		
		if (cid<0) {
			labs.push_back(C[i]);
			cid = labs.size();
		}
		C[i] = cid;		
	}
	for(int i = 0;i<N;i++){
		retC[i] = C[i];
		retX[i] = !!(X[i]);	
	}
	for(int i = 0;i < N;i++){
		adjList[i].clear();
		wList[i].clear();
	}
	adjList.clear();	
	wList.clear();	
	C.clear();	
	X.clear();	
}
