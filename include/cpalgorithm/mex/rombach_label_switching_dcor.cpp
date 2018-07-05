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
std::vector<int> ndord;
int N;
int M;
std::vector<int> deg; 
double alpha;
double beta;

double calc_coreness( int ord ){
	double c = 0.0;
	double bn = floor( N * beta);
	if(ord <= bn ){
		c = (1.0-alpha) / (2.0*bn) * (double)ord;
	}else{
		c = ((double)ord-bn)*(1.0-alpha)/(2*((double)N-bn))  + (1+alpha)/2.0;
	}	
	return c;
}


double rowSum_score( int nid, int sid){
	/*
	double retval = rowSumScore[nid];
	if( std::find(adjList[nid].begin(),adjList[nid].end(),sid) != adjList[nid].end()  ){
		retval-=calc_coreness( ndord[sid] );
	}
	return retval;
	*/
	double retval = 0;
	for( auto id : adjList[nid] ){
		if(id==sid) continue;
		retval+=calc_coreness( ndord[id] );
	}
	// Subtract rowSum_score on configuration model 
	for(int i = 0;i<N;i++){
		if(i==sid | i==nid) continue;
		retval-= (double)(deg[nid] * deg[i]) * calc_coreness( ndord[i] ) / (double) (2*M); 
	} 
	return retval;
}


double calc_swap_gain( int nid, int sid ){
	double c_nid = calc_coreness( ndord[nid] );
	double c_sid = calc_coreness( ndord[sid] );
	double rowSum_nid = rowSum_score(nid, sid);
	double rowSum_sid = rowSum_score(sid, nid);
	double dQ = (c_sid - c_nid) * rowSum_nid + (c_nid - c_sid) * rowSum_sid;
	dQ -=(c_nid*c_nid - c_sid*c_sid)*( deg[sid]*deg[sid]/(2*(double)M) )/2.0;
	dQ -=(c_sid*c_sid - c_nid*c_nid)*( deg[nid]*deg[nid]/(2*(double)M) )/2.0;
	return dQ; 
}
void swap_nodes( int nid, int nextnid ){
	// update ndord;	
	int tmp = ndord[nid];
	ndord[nid] = ndord[nextnid];
	ndord[nextnid]  = tmp;
	
}

std::vector<int> rombach_label_switching(){
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
	
	// ===========================
	// Initialise
	// ===========================
	std::shuffle(ndord.begin(), ndord.end(),mtrnd);
	
	// ===========================
	// Switching labels
	// ===========================
	std::vector<int> ord; ord.assign(N,0);
	for(int i = 0;i < N;i++) ord[i] = i;
		
	bool isupdated = false;
	int itmax = 100;int itnum=0;	
	do{
		isupdated = false;	
		std::shuffle(ord.begin(), ord.end(), mtrnd);
		for(int i = 0; i < N; i++){	
			int nid = ord[i]; // Get the id of node we shall update
			// calculate gain obtained by swapping label of nid and other labls
			int nextnid = nid; 
			double dQmax = -N;
			for( int sid =0;sid<N;sid++ ){
				if(sid==nid) continue;
				double dQ = calc_swap_gain(nid,sid);
				if(dQmax < dQ){
					nextnid = sid;
					dQmax = dQ;
				}
			}
			if( dQmax <= std::numeric_limits<double>::epsilon()) continue;
			isupdated = true;	
			swap_nodes(nid,nextnid);
		}
		itnum ++;
	}while( isupdated & itnum <itmax);
	//std::cout<<itnum<<std::endl;
	return ndord;
} 

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	N = (int) mxGetPr(prhs[1])[0];
   	M = (int) mxGetPr(prhs[2])[0]; //length of edge list
   	alpha = mxGetPr(prhs[3])[0];
   	beta = mxGetPr(prhs[4])[0];
	
	ndord.assign(N,0);
	for(int i = 0;i < N;i++) ndord[i] = i;
	
	for(int i = 0;i < N;i++){
		std::vector<int> tmp;
		adjList.push_back(tmp);
	}
	for(int i = 0;i<M;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+M]-1);
		adjList[rid].push_back(cid);
		adjList[cid].push_back(rid);
	}	
	deg.assign(N,0);
	for( int i = 0;i < N;i++ ){
		deg[i] = adjList[i].size();	
	} 
	rombach_label_switching();
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	double* R = mxGetPr(plhs[0]);
	for(int i = 0;i<N;i++){
		R[i] = calc_coreness( ndord[i] ) ;	
	}
	for(int i = 0;i < N;i++){
		adjList[i].clear();
	}
	adjList.clear();	
	ndord.clear();	
	deg.clear();	
}
