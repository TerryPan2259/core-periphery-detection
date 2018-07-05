#include "mex.h"
#include "be_cp.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
   	double M = (double) mxGetPr(prhs[2])[0]; //length of edge list
   	double* iscore =  mxGetPr(prhs[3]);
   	double pval = (double) mxGetPr(prhs[4])[0]; //length of edge list
	int numRandGraph = (int) mxGetPr(prhs[5])[0];
	
	vector<vector<int>> adjList;
	vector<bool> X;X.assign(N,false);
	for(int i = 0;i < N;i++){
		std::vector<int> tmp;
		adjList.push_back(tmp);
		if(iscore[i]>0.5){
			X[i]=true;
		}else{
			X[i]=false;
		}
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
	double estpval = 0;
	if(p>0.5){
		estpval =  estpval_boyd_stest_for_dense_net(adjList,X,N,M,p,numRandGraph,pval);
	}else{
		estpval =  estpval_boyd_stest(adjList,X,N,M,p,numRandGraph,pval);
	}

    	plhs[0] = mxCreateDoubleMatrix( (mwSize)1, (mwSize)1, mxREAL); 
    	double* R = mxGetPr(plhs[0]);
	R[0] = estpval;
}
