#include "mex.h"
#include "km_config_cont.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
   	int Enum = (double) mxGetPr(prhs[2])[0]; // length of edge list
	double* initC =  mxGetPr(prhs[3]);
	double* initX =  mxGetPr(prhs[4]);
		
	vector<vector<int>> adjList;
	vector<vector<double>> wList;
	vector<int> C;C.assign(N,0);
	vector<double> X;X.assign(N,0.0);
		
	for(int i = 0;i < N;i++){
		vector<int> tmp;
		adjList.push_back(tmp);
		vector<double> tmp2;
		wList.push_back(tmp2);
		C[i] = (int)(initC[i])-1;
		X[i] = initX[i];
	}
	
	double M = 0;
	for(int i = 0;i<Enum;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+(int)Enum]-1);
		double w = edgeList[i+2*(int)Enum];
		adjList[rid].push_back(cid);
		wList[rid].push_back(w);
		
		if(rid!=cid){
			adjList[cid].push_back(rid);
			wList[cid].push_back(w);
			M = M + w;
		}else{
			M = M + w;
		}
	}
	
	km_config_cont_label_switching(adjList, wList, C, X );
	
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
	
    	double* retC = mxGetPr(plhs[0]);

	for(int i = 0;i<N;i++){
		retC[i] = C[i]+1;
	}

}
