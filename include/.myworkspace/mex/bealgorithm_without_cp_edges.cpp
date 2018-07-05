#include "mex.h"
#include "be_cp.h"

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
	
	vector<bool> X;X.assign(N,false);
	bealgorithm_without_cp_edges(adjList, X, N, M, p, mtrnd );
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	double* R = mxGetPr(plhs[0]);
	for(int i = 0;i<N;i++) R[i] = !!(X[i]);
}
