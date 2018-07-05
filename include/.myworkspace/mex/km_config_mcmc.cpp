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

int weightedSampling(const vector<double>&S, mt19937_64& mtrnd){
	
	double sm = accumulate( S.begin(), S.end(), 0.0 );
	if(sm==0) return 0;
	uniform_real_distribution<double> udist(0.0,sm);
	double rv = udist(mtrnd);
	int idx = 0;double w = 0;
	for(auto it: S){
		w+=it;
		if( rv <= w ){
			return idx;
		}
		idx++;
	}	
}

double calcdQ( double dcore, double dperi, double degi, double Dcore, double Dperi, bool x, double M){
	return 2*( dcore + dperi*(!!(x)) - degi*(Dcore + Dperi*!!(x))/(2.0*M)) - !!(x)*degi*degi/(2.0*M); 
}

double calcQ(const vector<vector<int>>& adjList, const vector<vector<double>>& wList, const vector<double>& deg, const double M, const vector<int>& Ct, const vector<bool>& Xt){
	double retval = 0;
	int N = Ct.size();
	vector<double> Dc;Dc.assign(N,0.0);
	vector<double> Dp;Dp.assign(N,0.0);
	for(int i = 0;i < N;i++){
		for(int j = 0;j < adjList[i].size();j++){
			int nei = adjList[i][j];
			if( (Xt[i] | Xt[nei]) & (Ct[i] == Ct[nei] ) ) retval += wList[i][j];
		}
		Dc[Ct[i]]+=(double)!!(Xt[i])*deg[i];	
		Dp[Ct[i]]+=(double)!!(!Xt[i])*deg[i];	
	}
	
	for(int i =0;i<N;i++){
		retval-= Dc[i]*Dc[i]/(2*M) + 2*Dc[i]*Dp[i]/(2*M);	
	}
	return retval;	
}

//void bealgorithm( vector<int> indList, int N ){
void km_config_mcmc(const vector<vector<int>>& adjList, const vector<vector<double>>& wList, const double M, vector<int>& Cbest, vector<bool>& Xbest, int K, int maxNc, double beta, int maxStableLoopNum,int maxItNum){
	
	int N = adjList.size();
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
	vector<int> C;vector<bool> X;
	X.assign(N,true);C.assign(N,0);
	Xbest.assign(N,true);Cbest.assign(N,0);
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
	vector<double> deg;deg.assign(N,0.0);
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
	
		
	vector<double> toP;vector<double> toC;
	toC.assign(N,0.0);toP.assign(N,0.0);
	vector<double> Qcand;
	vector<int> Ccand;
	vector<bool> Xcand;
	vector<int> Cnei;
	bool isupdated = false;
	int itnum=0;int lastupdate = 1;
	double Q = calcQ(adjList, wList, deg, M, C, X);
	double Qbest = Q;	
	bool inCooling = false;
	while( itnum <=maxItNum ){
		if((itnum - lastupdate)>=maxStableLoopNum) {
			inCooling = true;
			beta*=10;
		}
			
		shuffle(ndord.begin(), ndord.end(),mtrnd);
		bool updateInCooling = false;
		for(int i = 0;i <N;i++){
			int nid = ndord[i];
			int ncid = C[nid];
			
			int neighbourNum = adjList[nid].size();
			fill(toC.begin(),toC.end(),0.0);	
			fill(toP.begin(),toP.end(),0.0);
			Cnei.clear();	
			int newcid = ncid;bool newx = X[nid];
			for(int j = 0;j<neighbourNum;j++){		
				int nei = adjList[nid][j];
				int cid = C[nei];
				toC[cid]+= (double)!!(X[nei]) * wList[nid][j];
				toP[cid]+= (double)!!(!X[nei]) * wList[nid][j];
				Cnei.push_back(cid);
			}
			sort( Cnei.begin(), Cnei.end() );
			Cnei.erase( unique( Cnei.begin(), Cnei.end() ), Cnei.end() );
			for(int j = 0;j<N;j++){
				if(find(Cnei.begin(),Cnei.end(),j)==Cnei.end()){
					Cnei.push_back(j);
					break;
				}
			}
		
				
			double dcore = Dcore[ncid] - deg[nid]*(double)!!(X[nid]);	
			double dperi = Dperi[ncid] - deg[nid]*(double)!!(!X[nid]);
			double dQold = calcdQ( toC[ncid], toP[ncid], deg[nid], dcore, dperi, X[nid], M );
			double Qcold = Qcc[ncid] - !!(X[nid])*(2*toC[ncid] -2*dcore * deg[nid]/(2*M) - deg[nid]*deg[nid]/(2*M));
			Qcand.clear();Ccand.clear();Xcand.clear();
			for(int j = 0;j<Cnei.size();j++){		
				int cid = Cnei[j];
					
				dcore = Dcore[cid] - deg[nid]*(double)!!(ncid==cid & X[nid]);	
				dperi = Dperi[cid] - deg[nid]*(double)!!(ncid==cid & !X[nid]);
		
				double dQc=calcdQ( toC[cid], toP[cid], deg[nid], dcore, dperi, true, M )-dQold;
				double dQp=calcdQ( toC[cid], toP[cid], deg[nid], dcore, dperi, false, M )-dQold;
				
				// check the density of the core subgraph 
				double dQcnew = 2*toC[cid] -2*dcore * deg[nid]/(2*M)-deg[nid]*deg[nid]/(2*M);
				bool movableToCore =  (Qcc[cid] + dQcnew >=0 ) & (Qcold >=0);
				bool movableToPeri =  (Qcc[cid]>=0)  & (Qcold >=0) | ( Ncore[cid]<=1);
				if(  movableToCore ){
					Qcand.push_back(exp(dQc*beta));
					Ccand.push_back(cid);
					Xcand.push_back(true);
				}
				if( movableToPeri ){
					Qcand.push_back(exp(dQp*beta));
					Ccand.push_back(cid);
					Xcand.push_back(false);
				}
			}
			int qid = weightedSampling(Qcand,mtrnd);
			newcid = Ccand[qid];
			newx = Xcand[qid];
			
			dcore = Dcore[newcid] - deg[nid]*(double)!!(ncid==newcid & X[nid]);	
			dperi = Dperi[newcid] - deg[nid]*(double)!!(ncid==newcid & !X[nid]);
			double dQ=calcdQ( toC[newcid], toP[newcid], deg[nid], dcore, dperi, newx, M )-dQold;
			double dQcore = !!(newx)*(2*toC[newcid] -2*dcore * deg[nid]/(2*M)-deg[nid]*deg[nid]/(2*M));
				
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
			C[nid] = newcid;	
			X[nid] = newx;	
			Q = Q+ dQ;	
			itnum++;	
			if(Q >Qbest){
				Qbest = Q;
				for(int i = 0;i<N;i++) {
					Cbest[i] = C[i];
					Xbest[i] = X[i];
				}
				if(inCooling) updateInCooling = true;
				lastupdate = itnum;
			}
		}
		if(updateInCooling==false & inCooling) break;
	}

} 



void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
   	int edgeNum = (int) mxGetPr(prhs[2])[0]; // length of edge list
   	int K = (int) mxGetPr(prhs[3])[0]; // maximum number of core-periphery pairs
   	int maxItNum = (int) mxGetPr(prhs[4])[0]; // maximum number of iterations
   	int maxStableLoopNum = (int) mxGetPr(prhs[5])[0]; // maximum number of intervals for updating Qmax 
	int maxNc = (int) mxGetPr(prhs[6])[0]; // maximum number of core nodes in a core-periphery pair
	double beta  = (double) mxGetPr(prhs[7])[0];
	
	
	vector<vector<int>> adjList;
	vector<vector<double>> wList;
	for(int i = 0;i < N;i++){
		std::vector<int> tmp;
		adjList.push_back(tmp);
		std::vector<double> tmp2;
		wList.push_back(tmp2);
	}
	double M = 0;
	for(int i = 0;i<edgeNum;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+edgeNum]-1);
		double w = edgeList[i+2*edgeNum];
		adjList[rid].push_back(cid);
		adjList[cid].push_back(rid);
		wList[rid].push_back(w);
		wList[cid].push_back(w);
		M = M + w;
	}
		
	
	
	vector<int> C;vector<bool> X;
	km_config_mcmc(adjList, wList, M, C, X, K, maxNc, beta, maxStableLoopNum,maxItNum);
	
	
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
		wList[i].clear();
	}
	adjList.clear();	
	wList.clear();	
	C.clear();	
	X.clear();	
}
