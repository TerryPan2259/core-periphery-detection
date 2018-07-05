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
std::vector<int> C;
std::vector<bool> X;
int N;
int M;
int K;
double p;

double calcQ(std::vector<int>Ct,std::vector<bool>Xt){
	double retval = 0;
	std::vector<double> Nc;Nc.assign(N,0.0);
	std::vector<double> Np;Np.assign(N,0.0);
	for(int i = 0;i < N;i++){
		for(int j = 0;j < adjList[i].size();j++){
			int nei = adjList[i][j];
			if( (Xt[i] | Xt[nei]) & (Ct[i] == Ct[nei] ) ) retval += wList[i][j];
		}
		Nc[Ct[i]]+=(double)!!(Xt[i]);	
		Np[Ct[i]]+=(double)!!(!Xt[i]);	
	}
	
	for(int i =0;i<N;i++){
		retval-= p*Nc[i]*(Nc[i]-1.0) + 2*p*Nc[i]*Np[i];	
	}
	return retval;	
}
//void bealgorithm( vector<int> indList, int N ){
void km_label_switching( ){
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
	std::vector<int> ndord;
	ndord.assign(N,0);
	for(int i = 0;i < N;i++) ndord[i] = i;
	
	C.assign(N,0);
	for(int i = 0;i < N;i++) {
		C[ ndord[i] ] = i % K ;
	}
		
	X.assign(N,true);
	
	std::vector<double> toP;std::vector<double> toC;
	
	std::uniform_real_distribution<double> udist(0.0,1.0);
	std::vector<int> Nc;std::vector<int> Np;
	toC.assign(N,0.0);toP.assign(N,0.0);
	Np.assign(N,0);Nc.assign(N,1);

	bool isupdated = false;
	do{
		isupdated = false;
		std::shuffle(ndord.begin(), ndord.end(),mtrnd);
	
		for(int i = 0;i <N;i++){
			int nid = ndord[i];
			int ncid = C[nid];
			
			int neighbourNum = adjList[nid].size();
			std::fill(toC.begin(),toC.end(),0.0);	
			std::fill(toP.begin(),toP.end(),0.0);
			
			int newcid = ncid;bool newx = X[nid];
			for(int j = 0;j<neighbourNum;j++){		
				int nei = adjList[nid][j];
				int w = wList[nid][j];
				int cid = C[nei];
				toC[cid]+= (double)!!(X[nei]) * w;
				toP[cid]+= (double)!!(!X[nei]) * w;
			}
			
			double dQ = 0;
			for(int j = 0;j<neighbourNum;j++){		
				int nei = adjList[nid][j];
				int cid = C[nei];
				//if( ncid==cid & X[nid] == X[nei] ) continue;
				double dQc = toC[cid] + toP[cid] - p * ( (double)Nc[cid] + (double)Np[cid] -  (double)!!(cid==ncid) );
				double dQp = toC[cid]  - p * ( (double)Nc[cid] - (double)!!(cid==ncid) );
			
				if( MAX(dQc,dQp)<dQ ) continue;	
				
				if (abs(dQc-dQp)<1e-8 ){ // if it is tie
					if( udist(mtrnd)<0.5){
						newx = true;
						newcid = cid;
						dQ = dQc;	
					}else{
						newx = false;
						newcid = cid;	
						dQ = dQp;	
					}
				}else{ // otherwise
					if( dQp < dQc ){
						newx = true;
						newcid = cid;
						dQ = dQc;	
					}else{
						newx = false;
						newcid = cid;	
						dQ = dQp;	
					}
				}
			}
			
			dQ = dQ - ( toC[ncid] + toP[ncid] -  p* ( (double)Nc[ncid] + (double)!!(X[nid])*(double)Np[ncid] -  (double)!!(X[nid]) ) );
				
			if(dQ<1e-6) continue;
			if( ncid==newcid & X[nid] == newx ) continue;
	
			if(X[nid]){
				Nc[ncid]--;
			}else{
				Np[ncid]--;
			}	
		
			if(newx){
				Nc[newcid]++;
			}else{
				Np[newcid]++;
			}
			isupdated = true;	
			C[nid] = newcid;	
			X[nid] = newx;	
		}
	}while( isupdated == true );

} 



void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	N = (int) mxGetPr(prhs[1])[0];
   	M = (int) mxGetPr(prhs[2])[0]; // length of edge list
   	K = (int) mxGetPr(prhs[3])[0]; // number of maximum core-periphery pair
	
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
	p = 2 * p / (N*(N-1));
		
	km_label_switching();
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
	}
	adjList.clear();	
	wList.clear();	
	C.clear();	
	X.clear();	
}
