#include <random>
#include <algorithm>
#include <vector>
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

double calcQ(
	const vector<vector<int>>& adjList, 
	const vector<vector<double>>& wList, 
	vector<bool>& S
	){
	int N = adjList.size();
	double Q = 0;
	for(int i = 0;i < N;i++){
		for(int j = 0;j < adjList[i].size();j++){
			if(S[i] != S[ adjList[i][j] ]){
				Q+=wList[i][j];
			}
		}
	}
	return Q;	
}

void group_bipartite_component(
	const vector<vector<int>>& adjList, 
	const vector<vector<double>>& wList, 
	vector<int>& C, 
	vector<bool>& S
	){
	
	int N = adjList.size();
	fill(C.begin(),C.end(),-1);
	vector<bool> x(N);fill(x.begin(),x.end(),false);
	vector<bool> xold(N);
	
	int cid = 0;
	int seednid = 0;	
	bool ischanged =false;
	do{
		x[seednid] = true;	
		xold[seednid] = true;	
		do{
			ischanged =false;
			
			for(int i = 0;i<N;i++){
				if(C[i]>=0) continue;	
				if(!xold[i]) continue;
				for(int j = 0;j<adjList[i].size();j++){
					int nei = adjList[i][j];
					
					if( ( S[i]==S[nei] )) continue; // ignore edges among the same side 
					
					if(xold[nei]==false){
						x[nei] = true;
						ischanged = true;	
					}
				}	
			}
			xold = x;	
		}while(ischanged);
		
		seednid = -1;	
		for(int i = 0;i<N;i++){
			if(x[i] & C[i]<0) C[i] = cid;
			if(C[i]<0 & seednid<0) seednid = i; //find the next seed node
		}
		cid++;
	}while(seednid>=0);
};

// flip 1 and 0 in S so that the number of ones are minimum 
void minimise_one(
	vector<int>& C, 
	vector<bool>& S
	){
	
	int N = C.size();
	int K = *max_element(C.begin(),C.end()) + 1;
	vector<int> L(K);fill(L.begin(),L.end(),0);
	vector<int> Ls(K);fill(Ls.begin(),Ls.end(),0);
	
	for(int i = 0;i<N;i++){
		L[ C[i] ]++;
		Ls[ C[i] ]+=!!(S[ i ]);
	}
	
	for(int i = 0;i<N;i++){
		if( L[C[i]] < 2 * Ls[C[i]] ) S[i] = !S[i]; // if one elements are greater than zeros
	}
		
};

void maxcut_mcmc(
	const vector<vector<int>>& adjList, 
	const vector<vector<double>>& wList, 
	vector<bool>& S, 
	double beta, 
	int maxStableLoopNum,
	int maxItNum
	){

	// --------------------------------
	// Initialise
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
	
	
	vector<bool> St=S;
	vector<int> ndord(N);
	double Q = calcQ(adjList,wList,S);
	for(int i = 0;i < N;i++) ndord[i] = i;
	double Qbest = Q;
		
	int itnum=0;int lastupdate = 1;
	while( itnum <=maxItNum & (itnum - lastupdate)<=maxStableLoopNum){
		shuffle(ndord.begin(), ndord.end(),mtrnd);
		for(int i = 0;i <N;i++){
			int nid = ndord[i];
				
			double toS = 0;double toR = 0;
			for(int j = 0;j<adjList[nid].size();j++){
				int nei = adjList[nid][j];
				toS+=(double)!!( St[ nei ] ) * wList[i][j];
				toR+=(double)!!( !St[ nei ] ) * wList[i][j];
			}
				
			double prob = exp( beta * toR )/( exp( beta * toS )   + exp( beta * toR ));
			
			if(udist(mtrnd)<=prob){
				Q-= toR*!!(St[nid]) + toS*!!(!St[nid]);
				Q+=toR;	
				St[nid] = true;	
			}else{
				Q-= toR*!!(St[nid]) + toS*!!(!St[nid]);
				Q+=toS;	
				St[nid] = false;	
			}	
			if(Q >Qbest){
				Qbest = Q;
				S = St;
				lastupdate = itnum;
			}
		}
		itnum++;	
	}
	Qbest = Q;
} 

void maxcut_kl(
	const vector<vector<int>>& adjList, 
	const vector<vector<double>>& wList, 
	vector<bool>& S
	){

	// --------------------------------
	// Initialise
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
		
	std::vector<bool>St(N);St = S;
	std::vector<double>toS(N);
	std::vector<double>toR(N);
	
	std::vector<bool>fixed(N);
	for( int j = 0;j < N;j++){
		fill(toS.begin(),toS.end(),0);fill(toR.begin(),toR.end(),0);	
		for(int i = 0;i< N;i++){
			for(int l = 0;l<adjList[i].size();l++){
				int nei = adjList[i][l];
				toS[i]+=(double)!!( St[ nei ] ) * wList[i][l];
				toR[i]+=(double)!!( !St[ nei ] ) * wList[i][l];
			}
		}
		
		std::fill(fixed.begin(),fixed.end(),false);
		double dQ = 0;
		double dQmax = -1 * std::numeric_limits<double>::max();
		int nid = 0;
		for( int i = 0;i < N;i++){
			double dqmax = -1 * std::numeric_limits<double>::max();
			
			// select a node of which we update the label 
			for(int k =0;k<N;k++){
				if( fixed[k] ) continue;
				double q = (toS[k]*!!(St[k]) + toR[k]*!!(!St[k])) - (toS[k]*!!(!St[k]) + toR[k]*!!(St[k]));
				if(dqmax < q ){
					nid = k;dqmax = q;
				}
			}

			
			St[ nid ] = !St[ nid ];
			for(int l = 0;l<adjList[nid].size();l++){
				int nei = adjList[nid][l];
				if(nid==nei) continue;
				toS[nei]+=(double)!!( !St[ nei ] ) * wList[nid][l] - (double)!!( St[ nei ] ) * wList[nid][l];
				toR[nei]+=(double)!!( St[ nei ] ) * wList[nid][l] - (double)!!( !St[ nei ] ) * wList[nid][l];
			}
			
			dQ = dQ + dqmax;
	
			//% Save the core-periphery pair if it attains the largest quality	
			if(dQmax < dQ){
				dQmax = dQ;
				if(dQmax>0) S = St;
			}
			//fixed( nid ) = true; % Fix the label of node nid
			fixed[ nid ] = true; //% Fix the label of node nid
		}
		 	
		if (dQmax <= std::numeric_limits<double>::epsilon()){
			break;
		}
		
		St = S;
	}
} 
