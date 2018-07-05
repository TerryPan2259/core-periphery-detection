#include "mex.h"
#include "km_config.h"
#include "km_config.cpp"

void init_random_number_generator(){
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
}


void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{
    double* edgeList = mxGetPr(prhs[0]);
    double* clist = mxGetPr(prhs[1]);
    double* xlist = mxGetPr(prhs[2]);
    int N = (int)mxGetPr(prhs[3])[0];
    int Enum = (int)mxGetPr(prhs[4])[0];
    int num_of_runs = (int)mxGetPr(prhs[5])[0];
    int num_of_rand_nets = (int)mxGetPr(prhs[6])[0];

    /* Parse Input */
    vector<vector<int> > A;
    vector<vector<double> > W;
    vector<int> c(N);
    vector<bool> x(N);
    for (int i = 0; i < N; i++) {
        vector<int> tmp;
        A.push_back(tmp);
        vector<double> tmp2;
        W.push_back(tmp2);
	c[i] = clist[i]-1;
	x[i] = xlist[i]==1;
    }

    for (int i = 0; i < Enum; i++) {
        int rid = (edgeList[i] - 1);
        int cid = (edgeList[i + (int)Enum] - 1);
        double w = edgeList[i + (int)2*Enum];
        A[rid].push_back(cid);
        W[rid].push_back(w);

        if (rid != cid) {
            W[cid].push_back(w);
            A[cid].push_back(rid);
        }
    }

    int K = *max_element(c.begin(),c.end()) + 1;
    vector<double> p_values;
    init_random_number_generator();
    cout<<"gen"<<endl;
    estimate_statistical_significance(A, W, c, x, num_of_runs, num_of_rand_nets, p_values);
    cout<<"end"<<endl;

    plhs[0] = mxCreateDoubleMatrix((mwSize)K, (mwSize)1, mxREAL);

    double* retPvals = mxGetPr(plhs[0]);

    for (int k = 0; k < K; k++) {
        retPvals[k] = p_values[k];
    }

}
