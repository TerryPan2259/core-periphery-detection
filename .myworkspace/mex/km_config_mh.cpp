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
using namespace std;


vector<vector<int>> adjList;
vector<vector<double>> wList;
vector<int> C;
vector<int> Cbest;
vector<bool> X;
vector<bool> Xbest;
vector<double> deg;
int N;
double M;
int K;
int maxNc;
double beta = 10;
int Type;

double calcQ(const vector<vector<int>>& adjList, const vector<vector<double>>& wList, const vector<double>& deg, const double M, const vector<int>& Ct, const vector<bool>& Xt){
	double retval = 0;
	int N = C.size();
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

double calcdQ( double dcore, double dperi, double degi, double Dcore, double Dperi, bool x ){
	if(Type==1){
		return 2*( dcore + dperi*(!!(x)) - degi*(Dcore + Dperi*!!(x))/(2.0*M)) - !!(x)*degi*degi/(2.0*M); 
	}else if(Type==2){
		return 2*( dcore + dperi*(2.0*!!(x)-1.0) - degi*(Dcore + Dperi*(2.0*!!(x)-1))/(2.0*M)) - (2*!!(x)-1)*degi*degi/(2.0*M); 
	}
}

double calcdQcore( double dcore, double degi, double Dcore, bool x ){
	return 2*( dcore - degi*Dcore/(2.0*M)) - degi*degi/(2.0*M); 
}

double calcQs( int cid ){
	double retval = 0;
	double Dcore = 0;
	for(int i = 0;i < N;i++){
		for(int j = 0;j < adjList[i].size();j++){
			int nei = adjList[i][j];
			if(C[i]==cid & C[nei]==cid & X[i] & X[nei]){ 
				retval+=wList[i][j];
			}	
		}	
		if(C[i]==cid & X[i]) Dcore+=deg[i];
	}
		
	retval = retval - (Dcore*Dcore)/(2*M);
	return retval;
}

//void bealgorithm( vector<int> indList, int N ){
void km_config_mcmc(int maxStableLoopNum,int maxItNum){
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
	vector<int> ndord;
	ndord.assign(N,0);
	for(int i = 0;i < N;i++) ndord[i] = i;
	
	C.assign(N,0);X.assign(N,true);
	vector<double> Nc;vector<double> Np;
	Np.assign(N,0.0);Nc.assign(N,0.0);
	for(int i = 0;i < N;i++) {
		C[ i ] = i;
		Nc[i]+=!!(X[i]);
		Np[i]+=!!(!X[i]);
		//C[ ndord[i] ] = i % K ;
		//Nc[i % K]+=!!(X[ndord[i]]);
		//Np[i % K]+=!!(!X[ndord[i]]);
		//if(Nc[i % K]+1>maxNc){
		//	Np[i % K]+=1.0;
		//	X[i]=false;
		//}else{
		//	Nc[i % K]+=1.0;
		//	X[i]=true;
		//}
	}
	
	vector<double> Dcore;Dcore.assign(N,0);
	vector<double> Dperi;Dperi.assign(N,0);
	vector<double> Qcc;Qcc.assign(N,0);
	deg.assign(N,0.0);
	for(int i = 0;i < N;i++) {
		for(int j = 0;j < wList[i].size();j++) {	
			deg[i]+=wList[i][j];
		}
	}
	for(int i = 0;i < N;i++) {
		Dcore[ C[i] ]+=(double)!!(X[i])*deg[i];
		Dperi[ C[i] ]+=(double)!!(!X[i])*deg[i];
		Qcc[ C[i] ]=-!!(X[i])*deg[i]*deg[i]/(2*M);
	}
		
	uniform_real_distribution<double> udist(0.0,1.0);
	uniform_int_distribution<int> idist(0,N-1);	
	
	double Q = calcQ(adjList, wList, deg, M, C, X);
	
	
	Cbest = C;Xbest = X;
	double Qbest = Q;
	int itnum=0;int lastupdate = 1;
	while( (itnum - lastupdate)<=maxStableLoopNum & itnum <=maxItNum ){
		// select a pair of nodes
		int nid = idist(mtrnd);
		int ncid = C[nid];
	
		// propose a pair id	
		int nei = nid;
		int idx = uniform_int_distribution<>(0,adjList[nid].size())(mtrnd);
		int newcid = -1;
		if(idx==adjList[nid].size()){ // propose a new core-periphery pair
			nei = nid;
			newcid = nid;
		}else{ // propose from neighbouring core-periphery pairs
			nei = adjList[nid][idx];
			newcid = C[nei];
		}
		
		// calculate the number of edges to the neighbouring core-periphery pairs 
		double toNewC = 0;double toNewP = 0;
		double toC = 0;double toP = 0;
		for(int j = 0;j<adjList[nid].size();j++){	
			int nei = adjList[nid][j];
			int cid = C[nei];
			if(newcid == cid){
				toNewC+=wList[nid][j]*!!(X[nei]);
				toNewP+=wList[nid][j]*!!(!X[nei]);	
			}
			if(ncid == cid){
				toC+=wList[nid][j]*!!(X[nei]);
				toP+=wList[nid][j]*!!(!X[nei]);	
			}
		}
		double Dcnew = Dcore[newcid] -  (double)(!!(newcid==ncid & X[nid]) )*deg[nid];
		double Dpnew = Dperi[newcid] -  (double)(!!(newcid==ncid & !X[nid]) )*deg[nid];
		double Dcold = Dcore[ncid] -  (double)(!!(X[nid]))*deg[nid];
		double Dpold = Dperi[ncid] -  (double)(!!(!X[nid]))*deg[nid];
		
		// propose a core or periphery assignment
		double dQcold = 2*toC -2*Dcold * deg[nid]/(2*M) -deg[nid]*deg[nid]/(2*M) ;	
		double dQcnew = 2*toNewC -2*Dcnew * deg[nid]/(2*M) -deg[nid]*deg[nid]/(2*M);
		double Qcold = Qcc[nid] - !!(X[nid])*dQcold;
			
		bool movableToCore =  (Qcc[newcid] + dQcnew >=0 ) & (Qcold >=0);
		bool movableToPeri =  (Qcc[newcid]>=0)  & (Qcold >=0) | ( Nc[newcid]<=N);
		
		bool newx = false;
		if( movableToCore & movableToPeri){
			newx = udist(mtrnd)<0.5;
		}else if( movableToPeri ){
			newx = false;
		}else{
			continue;
		}
	
		double Qcnew = Qcc[newcid] + !!(newx)*(dQcnew); 
		//cout<<nid<<" "<<adjList[nid].size()<<" "<<Nc[newcid]<<" "<<dQcore << tmpq <<" "<<newx<<endl;
			
		double dQ = calcdQ( toNewC, toNewP, deg[nid], Dcnew, Dpnew, newx )- calcdQ( toC, toP, deg[nid], Dcold, Dpold, X[nid] );
		
		if( udist(mtrnd) >= exp( beta*dQ/(2*M) ) ) continue;
		
		if(X[nid]){
			Dcore[ncid]-=deg[nid];
			Qcc[ncid]=Qcold;
			Nc[ncid]--;
		}else{
			Dperi[ncid]-=deg[nid];
			Qcc[ncid]=Qcold;
			Np[ncid]--;
		}
			
		if(newx){
			Dcore[newcid]+=deg[nid];
			Qcc[newcid]=Qcnew;
			Nc[newcid]++;
		}else{
			Dperi[newcid]+=deg[nid];
			Qcc[newcid]=Qcnew;
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
			cout<<Qbest<<endl;
		}
	}
} 



void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	N = (int) mxGetPr(prhs[1])[0];
   	int Enum = (double) mxGetPr(prhs[2])[0]; // length of edge list
   	K = (int) mxGetPr(prhs[3])[0]; // maximum number of core-periphery pairs
	K = N; //we fix
   	int maxItNum = (int) mxGetPr(prhs[4])[0]; // maximum number of iterations
   	int maxStableLoopNum = (int) mxGetPr(prhs[5])[0]; // maximum number of intervals for updating Qmax 
	maxNc = (int) mxGetPr(prhs[6])[0]; // maximum number of core nodes in a core-periphery pair
	beta  = (double) mxGetPr(prhs[7])[0];
	Type  = (int) mxGetPr(prhs[8])[0];
	for(int i = 0;i < N;i++){
		vector<int> tmp;
		adjList.push_back(tmp);
		vector<double> tmp2;
		wList.push_back(tmp2);
	}
	M =0;
	for(int i = 0;i<Enum;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+Enum]-1);
		double w = edgeList[i+2*Enum];
		adjList[rid].push_back(cid);
		adjList[cid].push_back(rid);
		wList[rid].push_back(w);
		wList[cid].push_back(w);
		M = M + w;
	}
		
	km_config_mcmc(maxStableLoopNum,maxItNum);
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	plhs[1] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
	
    	double* retC = mxGetPr(plhs[0]);
    	double* retX = mxGetPr(plhs[1]);
	
	vector<int> labs;
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
	for(int i = 0;i<N;i++){
		retC[i] = Cbest[i];
		retX[i] = !!(Xbest[i]);	
	}
	for(int i = 0;i < N;i++){
		adjList[i].clear();
		wList[i].clear();
	}
	adjList.clear();	
	wList.clear();	
	C.clear();	
	X.clear();	
	Cbest.clear();	
	Xbest.clear();	
	deg.clear();	
}
