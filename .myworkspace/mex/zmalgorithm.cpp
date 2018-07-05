#include <array>
#include <random>
#include <algorithm>
#include <bitset>
#include <vector>
#include <map>
#include <set>
#include "mex.h"
#include "math.h"
#include <random>
#include <iostream>     // std::cout, std::end
#include <stdlib.h>     /* abs */
#include <fstream>
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif


//void bealgorithm( vector<int> indList, int N ){
std::vector<bool> zmalgorithm( std::vector<int> indList, int N, int M, double maxItNum_bp, double maxItNum_em, double tol ){
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
		
	
	/* ============================	
	 Initialise 
	 ============================
	*/
	std::uniform_real_distribution<double> urand(0.0,1.0);
	std::vector<double> p; p.assign(4,0);
	std::vector<double> q; q.assign(N*2,0);
	std::vector<double> gamma; gamma.assign(2,0);
	std::map<int,std::vector<double> > eta, neweta; 
	p[0] = urand(mtrnd);p[2] = urand(mtrnd);
	p[1] = p[2];p[3] = urand(mtrnd);
	gamma[0] =urand(mtrnd);	gamma[1] =urand(mtrnd);
	double denom = gamma[0] + gamma[1];
	gamma[0]/=denom;gamma[1]/=denom;
	for(int i = 0; i < N ;i++){
		q[i] = urand(mtrnd); 
		q[i+N] = urand(mtrnd);
		denom = q[i] + q[i + N];
		q[i]/=denom;
		q[i+N]/=denom;
	}	
	for( auto ind : indList ){
		std::vector<double> etaij; etaij.assign(2,0);
		eta[ ind ] = etaij; 
		eta[ind][0] =  urand(mtrnd); eta[ind][1] =  urand(mtrnd);
		double denom = eta[ind][0] + eta[ind][1];
		eta[ind][0]/=denom;eta[ind][1]/=denom;
	}
	/* ============================	
	 EMalgorithm 
	 ============================
	*/
	int EMitnum = 0;int maxEMitnum = 100;
	double dif = 0;
	double dif2 = 0;
	std::vector<double> oldp; oldp = p; 
	std::vector<double> oldgamma; oldgamma = gamma;
	std::vector<double> newq; newq.assign(N*2,0);
	do{
		oldp = p;
		oldgamma = gamma;
		dif2 = 0;
		
		// ----------------------------------------------------------------
		// E-step
		// Beleaf propagation
		// ----------------------------------------------------------------
		// initialise eta	
		// update eta
		int itnum = 0;
		do{
			neweta = eta;
			dif = 0.0;	
			for( auto ind : indList ){
				int i = ind % N;
				int j = ind / N;
				for(int r = 0;r < 2;r++){
					double pd = 1;
					for(int k = 0;k<N;k++){	
						if (eta.find(i+k*N)!=eta.end()) continue;
						double qpsum = 0;
						for(int s = 0;s<2;s++){
							qpsum+=q[k + s*N] * p[ r + s*2 ];
						}
						pd = pd * MAX(1-qpsum, 0);
					}
					for(int k = 0;k<N;k++){	
						if (eta.find(i+k*N)==eta.end() | k==j) continue;
						double etapsum = 0;
						for(int s = 0;s<2;s++){
							etapsum+=eta[ k + N*i ][s] * p[r + s*2];
						}
						pd = pd * etapsum;
					}
					neweta[ind][r] = pd * gamma[r];
				}
				denom = neweta[ind][0] + neweta[ind][1];
				neweta[ind][0]/=denom;neweta[ind][1]/=denom;
				
				dif+=fabs(neweta[ind][0] - eta[ind][0]);
				dif+=fabs(neweta[ind][1] - eta[ind][1]);
			}
			// update q
			newq = q;
			for(int i = 0;i < N;i++){
				for(int r = 0;r < 2;r++){
					double pd = 1;
					for(int k = 0;k<N;k++){	
						if (eta.find(i+k*N)!=eta.end()) continue;
						double qpsum = 0;
						for(int s = 0;s<2;s++){
							qpsum+=q[k + s*N] * p[ r + s*2 ];
						}
						pd = pd * MAX(1-qpsum,0);
					}
					for(int k = 0;k<N;k++){	
						if (eta.find(i+k*N)==eta.end() ) continue;
						double etapsum = 0;
						for(int s = 0;s<2;s++){
							etapsum+=eta[ k + i*N ][s] * p[r + s*2];
						}
						pd = pd * etapsum;
					}
					newq[i + r*N] = gamma[r] * pd;	
				}
				double qir = newq[ i ];
				double qis = newq[ i + N ];
				newq[i] /= (qir + qis);
				newq[i + N] /= (qir + qis);
			
				dif+=fabs(newq[i] -q[i]);
				dif+=fabs(newq[i+N] -q[i +N ]);
			}	
				
			dif/=(double)( 2*N + 2*M );	
			
			eta = neweta;
			q = newq;
			itnum ++;
				
		}while( itnum < maxItNum_bp & dif > tol );
		

		// ----------------------------------------------------------------
		// M-step
		// ----------------------------------------------------------------
		dif2 = 0.0;
		gamma[0]= 0;gamma[1]= 0;
		for(int i = 0;i<N;i++){
			gamma[0]+= q[i];
			gamma[1]+= q[i + N ];
		}
		double denom = gamma[0] + gamma[1];
		gamma[0]/=denom;gamma[1]/=denom;
		dif2+=fabs(gamma[0] -oldgamma[0]);
		dif2+=fabs(gamma[1] -oldgamma[1]);
		
		std::vector<double> numer;numer.assign(4,0);
		std::vector<double> qrsij;qrsij.assign(4,0);
		for( auto ind : indList){
			int i = ind % N;
			int j = ind / N;
			
			int l = 0;denom = 0;
			for(int r = 0;r < 2;r++){
				for(int s = 0;s < 2;s++){
					qrsij[l]=eta[i+j*N][r] * eta[j+i*N][s] * p[r + s*2];
					denom+=qrsij[l];
					l++;
				}		
			}
			for(int ll = 0;ll < 4;ll++){
				numer[ll]+=qrsij[ll]/denom;
			}
		}
		std::vector<double> qrs; qrs.assign(2,0);
		for(int i = 0;i < N;i++){
			qrs[0]+=q[i];
			qrs[1]+=q[i + N];
		}	
		int l = 0;
		for(int r = 0;r < 2;r++){
			for(int s = 0;s < 2;s++){
				p[r + s*2] = MAX(MIN(1,numer[l] / ( qrs[r]  * qrs[s])),0);		
				dif2+=fabs(p[ r + 2*s] -oldp[r + 2*s]);
				l++;
			}
		}
		dif2/=(double)(2.0 + 4.0);
		EMitnum++;
	}while( dif2 > tol & EMitnum < maxItNum_em);
	//  --------------------------------------------- 
	std::vector<bool> C; C.assign(N,false);
	if( p[3] < p[0]  ){
		for(int i = 0;i < N;i++){
			if( q[i] > q[i + N]){
				C[i] = true;
			} 
		}	
	}else{
		for(int i = 0;i < N;i++){
			if( q[i] <  q[i + N]){
				C[i] = true;
			} 
		}	
	}
	return C;	
} 

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
   	int M = (int) mxGetPr(prhs[2])[0]; //length of edge list
   	double maxItNum_bp = mxGetPr(prhs[3])[0];
   	double maxItNum_em = mxGetPr(prhs[4])[0];
   	double tol = mxGetPr(prhs[5])[0];
	
	std::vector<int> L;
	for(int i = 0;i<M;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+M]-1);
		int ind = rid + cid * N;
		L.push_back(ind);
	}	
	std::vector<bool> C =zmalgorithm(L,N,M, maxItNum_bp,maxItNum_em,tol);
    	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)1, mxREAL); 
    	double* R = mxGetPr(plhs[0]);
	for(int i = 0;i<N;i++){
		if( C[i]){
			R[i] = 1;
		}else{
			R[i] = 0;
		}	
	}	
}
