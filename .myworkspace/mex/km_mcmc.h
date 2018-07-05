#include <array>
#include <algorithm>
#include <vector>
#include <random>
#include <iostream>     // std::cout, std::end
#include <fstream>
#include "math.h"
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif
using namespace std;

double calcQ(const vector<vector<int>>& adjList, const vector<vector<double>>& wList, const double p, int N, const vector<int>& Ct, const vector<bool>& Xt){
	double retval = 0;
	vector<double> Nc;Nc.assign(N,0.0);
	vector<double> Np;Np.assign(N,0.0);
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
void km_mcmc(const vector<vector<int>>& adjList, const vector<vector<double>>& wList, const double p, vector<int>& Cbest, vector<bool>& Xbest, int K, int maxNc, double beta, int maxStableLoopNum,int maxItNum){
	// --------------------------------
	// Initialise random generator 
	// --------------------------------
	int N = adjList.size();
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
	uniform_real_distribution<double> udist(0.0,1.0);
	uniform_int_distribution<int> idist(0,N-1);	
	
	// --------------------------------
	// Initialise arrays 
	// --------------------------------
	vector<int> C;
	vector<bool> X;
	vector<int> ndord;
	ndord.assign(N,0);
	for(int i = 0;i < N;i++) ndord[i] = i;
	
	C.assign(N,0);X.assign(N,true);
	vector<double> Nc;vector<double> Np;
	Np.assign(N,0.0);Nc.assign(N,0.0);
	for(int i = 0;i < N;i++) {
		C[ ndord[i] ] = i % K ;
		if(Nc[i % K]+1>maxNc){
			Np[i % K]+=1.0;
			X[i]=false;
		}else{
			Nc[i % K]+=1.0;
			X[i]=true;
		}
	}
	
		
	// --------------------------------
	// Initialise parameters 
	// --------------------------------
	double Q = 0.0;
	Q = calcQ(adjList, wList, p, N, C,X);
	Cbest = C;Xbest = X;
	
		
	// --------------------------------
	// MCMC 
	// --------------------------------
	double Qbest = Q;
	int itnum=0;int lastupdate = 1;
	while( (itnum - lastupdate)<=maxStableLoopNum & itnum <=maxItNum ){
		// select a pair of nodes
		int nid = idist(mtrnd);
		int ncid = C[nid];
	
		// propose a neighbouring core-perihpery pair or new pair	
		int nei = nid;
		int idx = uniform_int_distribution<>(0,adjList[nid].size())(mtrnd);
		int newcid = -1;
		if(idx==adjList[nid].size()){
			nei = nid;
			newcid = nid;
		}else{ // propose a neighbouring core-periphery pair
			nei = adjList[nid][idx];
			newcid = C[nei];
		}
			
		
		bool newx = X[nid];
		if(udist(mtrnd)<0.5  ){
			newx = true;
		}else{
			newx = false;
		}	
		//if (newcid == ncid ) newx = !X[nid];
		if (newx & Nc[newcid]+1>maxNc){
			newx = false;
		}
		
		// calculate the increment in Q	
		double dQ = 0.0;
		for(int j = 0;j<adjList[nid].size();j++){	
			int nei = adjList[nid][j];
			int cid = C[nei];
			dQ+=(double)(!!((newx | X[nei]) & newcid==cid))*wList[nid][j];
			dQ-=(double)(!!((X[nid] | X[nei]) & ncid==cid))*wList[nid][j];
		}
		double Ncnew = 0,Ncold=0,Npnew=0,Npold =0;
		Ncnew = Nc[newcid] -  (double)(!!(newcid==ncid & X[nid]) );
		Npnew = Np[newcid] -  (double)(!!(newcid==ncid & !X[nid]) );
		Ncold = Nc[ncid] -  (double)(!!(X[nid]) );
		Npold = Np[ncid] -  (double)(!!(!X[nid]) );
	
		dQ+= -p*( Ncnew + (double)(!!(newx))*Npnew);
		dQ+= +p*( Ncold + (double)(!!(X[nid]))*Npold);
		dQ = dQ*2.0;	
		
		if( udist(mtrnd) >= exp( beta*dQ )) continue;
		
			
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
		C[nid] = newcid;	
		X[nid] = newx;
		Q = Q + dQ;
		itnum++;	
		if(Q >Qbest){
			Qbest = Q;
			Cbest = C;
			Xbest = X;
			lastupdate = itnum;
		}
			
	}
	// remove redundant cp
	std::vector<int> labs;
	for(int i=0;i<N;i++){
		int cid = -1;
		for(int j=0;j<labs.size();j++){
			if(labs[j]==Cbest[i]){
				cid = j+1;
				break;
			}
		}
		
		if (cid<0) {
			labs.push_back(Cbest[i]);
			cid = labs.size();
		}
		Cbest[i] = cid;		
	}
} 



