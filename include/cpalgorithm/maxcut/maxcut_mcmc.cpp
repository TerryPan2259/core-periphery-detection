#include "mex.h"
#include "maxcut.h"

void mexFunction( int nlhs, mxArray *plhs[], 
	  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int Enum = (double) mxGetPr(prhs[1])[0]; // length of edge list
   	double* initS =  mxGetPr(prhs[2]);
   	int N = (int) mxGetPr(prhs[3])[0];
	int maxStableLoopNum = (int) mxGetPr(prhs[4])[0]; // Quality function Type 1: (c-c) + (c-p), Type 2: (c-c) + (c-p) - (p-p) 
	int maxItNum = (int) mxGetPr(prhs[5])[0]; // Quality function Type 1: (c-c) + (c-p), Type 2: (c-c) + (c-p) - (p-p) 
	double beta = (double) mxGetPr(prhs[6])[0]; // Quality function Type 1: (c-c) + (c-p), Type 2: (c-c) + (c-p) - (p-p) 
	
	
	vector<vector<int>> adjList;
	vector<vector<double>> wList;
	vector<bool> S(N);
	vector<int> C(N);
		
	for(int i = 0;i < N;i++){
		vector<int> tmp;
		adjList.push_back(tmp);
		vector<double> tmp2;
		wList.push_back(tmp2);
		S[i] = ((int)initS[i])==1;
	}
	
	for(int i = 0;i<Enum;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+(int)Enum]-1);
		double w = edgeList[i+2*(int)Enum];
		adjList[rid].push_back(cid);
		wList[rid].push_back(w);
		adjList[cid].push_back(rid);
		wList[cid].push_back(w);
	}
	
		
	maxcut_mcmc( adjList, wList, S, beta, maxStableLoopNum, maxItNum );
	group_bipartite_component( adjList, wList, C, S );
	minimise_one(C,S);

	vector<vector<double>> cplist;	
	int cid = 0;
	for(int i = 0;i<N;i++){
		if(S[i]){
			for(int j = 0;j<adjList[i].size();j++){
				int nei = adjList[i][j];
				if(S[nei]) continue;
				
				vector<double>tmp(4);
				tmp[0] = i; // core 	
				tmp[1] = nei; // peri.
				tmp[2] = cid; // cpid
				tmp[3] = wList[i][j]; // strength	
				cplist.push_back(tmp);
			}
			cid++;
		}
	}

	int L = cplist.size();	
	
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)L, (mwSize)4, mxREAL); 
    	plhs[1] = mxCreateDoubleMatrix( (mwSize)1, (mwSize)1, mxREAL); 
	
    	double* retS = mxGetPr(plhs[0]);
	double* Q  = mxGetPr(plhs[1]);
	Q[0] = calcQ(adjList,wList,S);
	
	for(int i = 0;i<L;i++){
		retS[i] = cplist[i][0]+1;	
		retS[i+L] = cplist[i][1]+1;	
		retS[i+2*L] = cplist[i][2]+1;	
		retS[i+3*L] = cplist[i][3];	
	}

}
