#include <array>
#include <random>
#include <algorithm>
#include <bitset>
#include <vector>
#include <set>
#include "mex.h"
#include "math.h"
#include <random>
#include <iostream>     // cout, end
#include <fstream>
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif
using namespace std;

vector<vector<int>> adjList;
vector<vector<double>> wList;
vector<int> C;
vector<bool> X;
int N;
double M;
int K;
int Type;
vector<double> deg;
double calcdQ( double dcore, double dperi, double degi, double Dcore, double Dperi, bool x ){
	if(Type==1){
		return 2*( dcore + dperi*(!!(x)) - degi*(Dcore + Dperi*!!(x))/(2.0*M)) - !!(x)*degi*degi/(2.0*M); 
	}else if(Type==2){
		return 2*( dcore + dperi*(2.0*!!(x)-1.0) - degi*(Dcore + Dperi*(2.0*!!(x)-1))/(2.0*M)) - (2*!!(x)-1)*degi*degi/(2.0*M); 
	}
}

double calcdQcore( double dcore, double degi, double Dcore ){
	return 2*( dcore - degi*(Dcore)/(2.0*M)) - degi*degi/(2.0*M); 
}

//void bealgorithm( vector<int> indList, int N ){
void km_config_label_switching(){
	
	//% --------------------------------
	//% Initialise
	//% --------------------------------
	mt19937_64 mtrnd;
    	int seeds[624];
       	size_t size = 624*4; //Declare size of data
       	ifstream urandom("/dev/urandom", ios::in | ios::binary); //Open stream
       	if (urandom) //Check if stream is open
       	{
       	    urandom.read(reinterpret_cast<char*>(seeds), size); //Read from urandom
       	    urandom.close(); //close stream
       	}
       	else //Open failed
       	{
            		cerr << "Failed to open /dev/urandom" << endl;
       	}
    	seed_seq seed(&seeds[0], &seeds[624]);
    	mtrnd.seed(seed);
	
	
	// ===========================
	// Initialise
	// ===========================
	X.assign(N,true);
	C.assign(N,0);
	vector<int> ndord;
	ndord.assign(N,0);
	for(int i = 0;i < N;i++) ndord[i] = i;
	for(int i = 0;i < N;i++) {
		C[ i ] = i ;
	}
	
	vector<double> Dcore;Dcore.assign(N,0);
	vector<double> Dperi;Dperi.assign(N,0);
	vector<int> Ncore;Ncore.assign(N,0);
	vector<int> Nperi;Nperi.assign(N,0);
	deg.assign(N,0.0);
	vector<double> Qcc;Qcc.assign(N,0.0); // contribution of core-core connections to Q
	for(int i = 0;i < N;i++) {
		for(int j = 0;j < wList[i].size();j++) {	
			deg[i]+=wList[i][j];
		}
	}
	for(int i = 0;i < N;i++) {
		Dcore[ C[i] ]+=(double)!!(X[i])*deg[i];
		Dperi[ C[i] ]+=(1.0-(double)!!(X[i]))*deg[i];
		Ncore[ C[i] ]+=(double)!!(X[i]);
		Nperi[ C[i] ]+=(1.0-(double)!!(X[i]));
		Qcc[C[i]] = -deg[i]*deg[i]/(2*M)*!!(X[i]);
	}
	
		
	uniform_real_distribution<double> udist(0.0,1.0);
	vector<double> toP;vector<double> toC;
	toC.assign(N,0.0);toP.assign(N,0.0);
	bool isupdated = false;
	int itNum = 0;
	bool inCooling = false;
	do{
		isupdated = false;
		shuffle(ndord.begin(), ndord.end(),mtrnd);
		for(int i = 0;i <N;i++){
			int nid = ndord[i];
			int ncid = C[nid];
			
			int neighbourNum = adjList[nid].size();
			fill(toC.begin(),toC.end(),0.0);	
			fill(toP.begin(),toP.end(),0.0);
				
			int newcid = ncid;bool newx = X[nid];
			for(int j = 0;j<neighbourNum;j++){		
				int nei = adjList[nid][j];
				int cid = C[nei];
				toC[cid]+= (double)!!(X[nei]) * wList[nid][j];
				toP[cid]+= (double)!!(!X[nei]) * wList[nid][j];
			}
			
			double dcore = Dcore[ncid] - deg[nid]*(double)!!(X[nid]);	
			double dperi = Dperi[ncid] - deg[nid]*(double)!!(!X[nid]);
			double dQold = calcdQ( toC[ncid], toP[ncid], deg[nid], dcore, dperi, X[nid] );
			double dQ = 0;
			double dQcore = 0;
			double dQcold = 2*toC[ncid] -2*dcore * deg[nid]/(2*M) - deg[nid]*deg[nid]/(2*M);	
			double Qcold = Qcc[ncid] - !!(X[nid])*dQcold;
			for(int j = 0;j<neighbourNum;j++){		
				int nei = adjList[nid][j];
				int cid = C[nei];
					
				dcore = Dcore[cid] - deg[nid]*(double)!!(ncid==cid & X[nid]);	
				dperi = Dperi[cid] - deg[nid]*(double)!!(ncid==cid & !X[nid]);
		
				double dQc=calcdQ( toC[cid], toP[cid], deg[nid], dcore, dperi, true )-dQold;
				double dQp=calcdQ( toC[cid], toP[cid], deg[nid], dcore, dperi, false )-dQold;
			
				if( MAX(dQc,dQp)<dQ ) continue;
				if( MAX(dQc,dQp)<0 ) continue;
					
				// check the density of the core subgraph 
				double dQcnew = 2*toC[cid] -2*dcore * deg[nid]/(2*M)-deg[nid]*deg[nid]/(2*M);
				bool movableToCore =  (Qcc[cid] + dQcnew >=0 ) & (Qcold >=0);
				bool movableToPeri =  (Qcc[cid]>=0)  & (Qcold >=0) | ( Ncore[cid]<=1);
				
				if(abs(dQp-dQc)<1e-8){
					newx = udist(mtrnd)<0.5;
					if( newx & movableToCore ){
						newx = true;
						newcid = cid;
						dQ = dQc;	
						dQcore = dQcnew;
					}else if(newx==false & movableToPeri){
						newx = false;
						newcid = cid;	
						dQ = dQp;	
						dQcore = 0;
					}
				}else{
					if( dQp < dQc & movableToCore ){
						newx = true;
						newcid = cid;
						dQ = dQc;	
						dQcore = dQcnew;
					}
					if( dQp > dQc & movableToPeri ){
						newx = false;
						newcid = cid;	
						dQ = dQp;	
						dQcore = 0;
					}
				}
			}
			
			if(dQ<1e-9) continue;
				
			if( ncid==newcid & X[nid] == newx ) continue;
				
			if(X[nid]){
				Dcore[ncid]-=deg[nid];
				Ncore[ncid]--;
				Qcc[ncid]=Qcold;
			}else{
				Qcc[ncid]=Qcold;
				Dperi[ncid]-=deg[nid];
				Nperi[ncid]--;
			}
				
			if(newx){
				Dcore[newcid]+=deg[nid];
				Qcc[newcid]+=dQcore;
				Ncore[newcid]++;
			}else{
				Dperi[newcid]+=deg[nid];
				Nperi[newcid]++;
				Qcc[newcid]+=dQcore;
			}
			isupdated = true;	
			C[nid] = newcid;	
			X[nid] = newx;	
		}
		itNum++;
	}while( isupdated == true );

} 



void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	N = (int) mxGetPr(prhs[1])[0];
   	int Enum = (double) mxGetPr(prhs[2])[0]; // length of edge list
   	K = (int) mxGetPr(prhs[3])[0]; // number of maximum core-periphery pair
   	//K = N; // we fix 
	Type = (int) mxGetPr(prhs[4])[0]; // Quality function Type 1: (c-c) + (c-p), Type 2: (c-c) + (c-p) - (p-p) 
		
	for(int i = 0;i < N;i++){
		vector<int> tmp;
		adjList.push_back(tmp);
		vector<double> tmp2;
		wList.push_back(tmp2);
	}
	M = 0;
	for(int i = 0;i<Enum;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+(int)Enum]-1);
		double w = edgeList[i+2*(int)Enum];
		adjList[rid].push_back(cid);
		adjList[cid].push_back(rid);
		wList[rid].push_back(w);
		wList[cid].push_back(w);
		M = M + w;
	}
	km_config_label_switching();
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	plhs[1] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
	
    	double* retC = mxGetPr(plhs[0]);
    	double* retX = mxGetPr(plhs[1]);
	
	vector<int> labs;
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
	deg.clear();	
}
