#include "mex.h"
#include "bp.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
   	int Enum = (double) mxGetPr(prhs[2])[0]; // length of edge list
		
	vector<vector<int>> eList;
	vector<int> C;C.assign(N,0.0);
	double M = 0;
	for(int i = 0;i<Enum;i++){
		vector<int>edge(2);
		edge[0] = (int) edgeList[i]-1; 
		edge[1] = (int) edgeList[i+(int)Enum]-1; 
		eList.push_back(edge);
	}
	
	bp(eList, C, N, 2, 10, 0.000001);
	
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
	
    	double* retC = mxGetPr(plhs[0]);

	for(int i = 0;i<N;i++){
		retC[i] = C[i]+1;
	}

}
