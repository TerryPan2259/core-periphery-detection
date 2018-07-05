#include <stdio.h>
#include <random>
#include <stdlib.h> 
#include <string.h>
#include <iostream>
#include <math.h> 
#include <fstream>
#include <limits>
#include <ctime>
#include <map>
#include <vector>
#include <algorithm>
//#include <gsl/gsl_rng.h>  //random number generator for initialization.

#if !defined(MAX_K)
#define	MAX_K 100
#endif

using namespace std;

//Constants for the graph.
//const int N = 34; //number of nodes
//const int M = 78;//number edges
//const int K = 2; //number of groups
//const int max_ite = 50; //number of maximum iteration
//const long double tol = 0.000001; //tolerance, not implemented in this version, could be defined by the user.

//Variable declaration.
//long int *adjList[N]; //adjacency list of the network.
//long int Degree[N];
//long int EdgeList[M][2]; //edgelist

//long double Gamma[K];
//long double Omega[K][K];
//long double One_point[N][K]; //one point marginal

//map<pair<int,int>, array<array<long double,2>, 2> > Two_point; //two point marginal
//map<pair<int,int>, array<long double,2> > Cavity; //message matrix, stored the edges as keys

//double External_field[K]; //external field term as sum of logs


//Function declaration.
void update_Gamma(vector<long double>& Gamma, const int N, const int K, vector<vector<long double>>& One_point);
void update_Omega(vector<vector<long double>>& Omega,  const int N, const int K,const vector<vector<int>>& adjList,  const vector<int>& Degree, const vector<long double>& Gamma, map<pair<int,int>, array<array<long double,MAX_K>, MAX_K> >& Two_point);

void initilize_params(const vector<vector<int>>& edgeList, const int N, const int K, vector<vector<int>>& adjList, vector<int>& Degree, vector<long double>& Gamma, vector<vector<long double>>& Omega, vector<vector<long double>>& One_point, map<pair<int,int>, array<array<long double,MAX_K>, MAX_K> >& Two_point, map<pair<int,int>, array<long double,MAX_K> >& Cavity, vector<double>& External_field);

void compute_external(vector<double>& External_field, const int N, const int K, const vector<vector<long double>>& Omega,  const vector<vector<long double>>& One_point);
void predict(vector<int>&Glist, const int N, const int K, vector<vector<long double>>& One_point);
void compute_one_point(vector<vector<long double>>& One_point,  int N,  int K,  vector<vector<int>>& adjList,  vector<int>& Degree,  vector<long double>& Gamma,  vector<vector<long double>>& Omega,  map<pair<int,int>, array<long double,MAX_K> >& Cavity,  vector<double>& External_field);
void compute_two_point(map<pair<int,int>, array<array<long double,MAX_K>, MAX_K> >& Two_point, const int N, const int K, const vector<vector<int>>& adjList, const vector<int>& Degree, const vector<vector<long double>>& Omega,  map<pair<int,int>, array<long double,MAX_K> >& Cavity);
void compute_cavity(map<pair<int,int>, array<long double,MAX_K> >& Cavity, const int N, const int K, const vector<vector<int>>& adjList, const vector<int>& Degree, const vector<long double>& Gamma, const vector<vector<long double>>& Omega, const vector<vector<long double>>& One_point, const vector<double>& External_field);

// modified
void bp(const vector<vector<int>>& edgeList, vector<int>& C, const int N, const int K, const int max_ite, const long double tol);


void bp(const vector<vector<int>>& edgeList, vector<int>& C, const int N, const int K, const int max_ite, const long double tol)
{

    //initialization
    vector<vector<int>> adjList;
    vector<int> Degree;
    vector<long double> Gamma;
    vector<vector<long double>> Omega;
    vector<vector<long double>> One_point;
    map<pair<int,int>, array<array<long double,MAX_K>, MAX_K> > Two_point; //two point marginal
    map<pair<int,int>, array<long double,MAX_K> > Cavity; //message matrix, stored the edges as keys
    vector<double> External_field; //external field term as sum of logs
    
    initilize_params(edgeList, N, K, adjList, Degree, Gamma, Omega, One_point, Two_point, Cavity, External_field);

    //main loop
    int ite = 0;
    while (ite < max_ite)
    {
        //update message and two point marginal
        compute_external(External_field, N, K, Omega, One_point);                
        compute_one_point(One_point, N, K, adjList, Degree, Gamma, Omega, Cavity, External_field);      
        compute_cavity(Cavity, N, K, adjList, Degree, Gamma, Omega, One_point, External_field);
	
        //update parameters
        compute_two_point(Two_point, N, K, adjList, Degree, Omega,  Cavity);
        update_Omega(Omega, N, K, adjList, Degree, Gamma, Two_point);
        update_Gamma(Gamma, N, K, One_point);    

        ite++;
    }

    fill(C.begin(),C.end(),0);
    predict(C, N, K, One_point);
}

//based on the implementation by Brian Karrer.

//update the parameters based on current belief
void update_Gamma(vector<long double>& Gamma, const int N, const int K, vector<vector<long double>>& One_point)
{
    int r;
    long int i;
    long double NORM=0;    

    fill(Gamma.begin(), Gamma.end(), 0);
    for(r = 0; r < K; r++)
    {
        for(i = 0; i < N; i++)
        {
            Gamma[r] += One_point[i][r]; 
        }
        NORM += Gamma[r];        
    }
    
    //normalize
    for(r=0; r < K; r++)
    {
        Gamma[r] = Gamma[r]/NORM;
    }   
/* 
    cout << "The Gamma is now" << endl;
    cout << Gamma[0] << " " << Gamma[1] << endl;
*/
	return;
}

void update_Omega(vector<vector<long double>>& Omega,  const int N, const int K, const vector<vector<int>>& adjList,const vector<int>& Degree, const vector<long double>& Gamma, map<pair<int,int>, array<array<long double,MAX_K>, MAX_K> >& Two_point )
{
    long int i,j;
    int r,s,it;
    long double SIGMA;

    for(r = 0; r < K; r++)
    {
        for(s = 0; s < K; s++)
        {
            SIGMA = 0;
            for(i = 0; i < N; i++)
            {
                for(it = 0; it < Degree[i]; it++)
                {
                    j = adjList[i][it];
                    SIGMA += Two_point[make_pair(i,j)][r][s];
                }
            }
            //update the mixing matrix.
            Omega[r][s] = SIGMA/(N*N*Gamma[r]*Gamma[s]);
        }
    } 
/*
    cout << "The mixing matrix is now" << endl;
    cout << Omega[0][0] << " " << Omega[0][1] << endl;
    cout << Omega[1][0] << " " << Omega[1][1] << endl;
*/
	return;
}

//make prediction based on belief. i,e. assign node to the group with highest belief.
void predict(vector<int>&Glist, const int N, const int K, vector<vector<long double>>& One_point)
{   
    long int i;
    for (i = 0; i < N; i++)
    {
        Glist[i] = distance(One_point[i].begin(),max_element(One_point[i].begin(),One_point[i].end()));        
        //Glist[i] = distance(One_point[i],max_element(One_point[i].begin(),One_point[i].end()));        
    }    
	return;
}

void initilize_params(const vector<vector<int>>& edgeList, const int N, const int K, vector<vector<int>>& adjList, vector<int>& Degree, vector<long double>& Gamma, vector<vector<long double>>& Omega, vector<vector<long double>>& One_point, map<pair<int,int>, array<array<long double,MAX_K>, MAX_K> >& Two_point, map<pair<int,int>, array<long double,MAX_K> >& Cavity, vector<double>& External_field)
{
     Degree.assign(N,0);
     for(int i=0; i < N; i++)
     {	
	vector<int> tmp;
	adjList.push_back(tmp);
			
	vector<long double> tmp2(K);
	One_point.push_back(tmp2);
     }
     for(int i=0; i < K; i++)
     {	
	vector<long double> tmp(K);
	Omega.push_back(tmp);
     }
     Gamma.assign(K,0.0);
     External_field.assign(K,0.0);
     
     // First we count the degrees by scanning through the list once
     for(int i=0; i < edgeList.size(); i++)
     {
	int src = edgeList[i][0];
	int dst = edgeList[i][1];
	adjList[src].push_back(dst);
	adjList[dst].push_back(src);
     }

     double avedeg = 0;	     
     for(int i=0; i < N; i++)
     {	
	Degree[i] = adjList[i].size();
	avedeg+=(double)Degree[i];
     }
	avedeg/=(double)N;

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

	//The initialization can be achieved many ways and different initialization 
	//will result in possibly different results. Try different initialization 
	//point as well as initialization method for better results. Here the 
	//method is the 'apriori' method where we have a rough idea of what the parameters
	//should be.
	double p = (avedeg * (double) N)/ (double)(N*(N-1));
	Gamma[0] = 0.5;
	Gamma[1] = 0.5;
	Omega[0][0] = 0.5;
	Omega[0][1] = p*p;
	Omega[1][0] = p*p;
	Omega[1][1] = p*p*p*p;
	
	//one could randomly generate a,b repeatedly for each node/edge. one value usually work just as well.
	srand(time(NULL)); //random seed 
	long double a;		
	long double b;
	a = ((long double) rand() / (RAND_MAX));
	b = ((long double) rand() / (RAND_MAX));	

    //initialize one point marginal to unbalanced belief.
    for(int i = 0; i < N ;i++)
    {
       One_point[i][0] = !!(Degree[i]>=avedeg);
       One_point[i][1] = 1-!!(Degree[i]>=avedeg);
    }	    
/*
	for(int i = 0;i < K;i++){
		Gamma[i] = 1.0/(double)(K);
		for(int j = 0;j < K;j++){
			Omega[i][j] = udist(mtrnd);
		}
	}
*/
    //initialize the message on the edges; also initialize the two point marginal to zeros.
    for(int i = 0; i < edgeList.size(); i++)
    {
        for(int r = 0; r < K; r++)
        {            
            Cavity[make_pair(edgeList[i][0],edgeList[i][1])][r] = r==0 ? b:1-b;
            Cavity[make_pair(edgeList[i][1],edgeList[i][0])][r] = r==0 ? b:1-b;
            Two_point[make_pair(edgeList[i][1],edgeList[i][0])] = {{{0,0},{0,0}}};
            Two_point[make_pair(edgeList[i][0],edgeList[i][1])] = {{{0,0},{0,0}}};            
        }        
    }    
}

//compute the external field term 
void compute_external(vector<double>& External_field, const int N, const int K, const vector<vector<long double>>& Omega,  const vector<vector<long double>>& One_point )
{
    long int i;
    int r,s;
    long double SIGMA;

    fill(External_field.begin(),External_field.end(),0);
    
    for(r = 0; r < K; r++)
    {
        for(i=0; i < N; i++)
        {
            SIGMA = 0;
            for(s = 0; s < K; s++)
            {
                SIGMA += One_point[i][s]*exp(-Omega[r][s]);
            }       
            External_field[r] += log(SIGMA);                  
        }
    }           
    return;      
}

//compute the one point marginal
void compute_one_point(vector<vector<long double>>& One_point,  int N,  int K,  vector<vector<int>>& adjList,  vector<int>& Degree,  vector<long double>& Gamma,  vector<vector<long double>>& Omega,  map<pair<int,int>, array<long double,MAX_K> >& Cavity,  vector<double>& External_field)
{
    long double BUFFER[K];
    long int i,j;
    int r,s,it;
    long double SIGMA1,SIGMA2,PROD;
    
    for(i = 0; i < N; i++)
    {
        fill(BUFFER,BUFFER+K,0);
        for(r = 0; r < K; r++)
        {
            PROD = 0;
            for(it = 0; it < Degree[i]; it++)
            {
                j = adjList[i][it];
                SIGMA1 = 0;
                SIGMA2 = 0;
                for(s = 0; s < K; s++)
                {                       
                    SIGMA1 += Cavity[make_pair(j,i)][s]*Omega[r][s]*exp(-Omega[r][s]);
                    SIGMA2 += One_point[j][s]*exp(-Omega[r][s]);                    
                }                
                PROD += log(SIGMA1)- log(SIGMA2);                    
            }
            BUFFER[r] = log(Gamma[r]) + PROD + External_field[r];
        }
        //normalize and set the value
        long double x = 1/(exp(BUFFER[0]-BUFFER[1])+1);        
        One_point[i][0] = 1 -x;
        One_point[i][1] = x;      
/*
            //normalize
            long double denom = 0;
	    long double md = 0; 
            for(int l =0;l<K;l++){
		md +=BUFFER[l];
	    }
	    md/=(long double) K;
            for(int l =0;l<K;l++){
		denom+= exp(BUFFER[l]-md);
	    }
            for(int l =0;l<K;l++){
                One_point[i][l] = exp(BUFFER[l]-md) / denom;
	    }
*/ 
    }
        
    return;
}

//compute the two point marginal
void compute_two_point(map<pair<int,int>, array<array<long double,MAX_K>, MAX_K> >& Two_point, const int N, const int K, const vector<vector<int>>& adjList, const vector<int>& Degree, const vector<vector<long double>>& Omega,  map<pair<int,int>, array<long double,MAX_K> >& Cavity)
{       
    long double BUFFER[K][K];
    long double NORM,PROD; //normalization term, sum of the individual q_{ij}^{rs}.
    long int i,j;
    int r,s,it;

    for(i = 0; i < N; i++)
    {
        for(it = 0; it < Degree[i]; it++)
        {
            j = adjList[i][it];
            NORM = 0;
            fill(BUFFER[0],BUFFER[0]+K,0);
            fill(BUFFER[1],BUFFER[1]+K,0);
            //summing over the four possible values of two point marginal
            for(r = 0; r < K; r++)
            {
                for(s = 0; s < K; s++)
                {
                    PROD = Omega[r][s]*exp(Omega[r][s])*Cavity[make_pair(i,j)][r]*Cavity[make_pair(j,i)][s];
                    NORM += PROD;
                    BUFFER[r][s] = PROD;                    
                }
            }
            
            //normalizing the two point marginal so they sum to unity
            for(r = 0; r < K; r++)
            {
                for(s = 0; s < K; s++)
                {
                    Two_point[make_pair(i,j)][r][s] = BUFFER[r][s]/NORM;
                }
            }                           
        }
    }
    return;
}

//compute the message matrix
//only update the edge terms, the non-egde terms will just be the value of the one point marginal
void compute_cavity(map<pair<int,int>, array<long double,MAX_K> >& Cavity, const int N, const int K, const vector<vector<int>>& adjList, const vector<int>& Degree, const vector<long double>& Gamma, const vector<vector<long double>>& Omega, const vector<vector<long double>>& One_point, const vector<double>& External_field)
{
    long double BUFFER[K];
    long int i,j,k;
    int r,s,it1,it2; //use two iterators for two edges
    long double SIGMA1,SIGMA2,PROD;

    for(i = 0; i < N; i++)
    {
        for(it1 = 0; it1 < Degree[i]; it1++)
        {
            j = adjList[i][it1];
            for(r = 0; r < K; r++)
            {
                PROD = 0;
                for(it2 = 0; it2 < Degree[i]; it2++)
                {
                    k = adjList[i][it2];
                    if(k != j)
                    {   
                        SIGMA1 = 0;
                        SIGMA2 = 0;
                        for(s = 0; s < K; s++)
                        {
                            SIGMA1 += Cavity[make_pair(k,i)][s]*Omega[r][s]*exp(-Omega[r][s]);
                            SIGMA2 += One_point[k][s]*exp(-Omega[r][s]);
                        }
                        PROD += log(SIGMA1) - log(SIGMA2);
                    }                                        
                }
                BUFFER[r] = log(Gamma[r]) + PROD + External_field[r];
            }
            //normalize
            long double x = 1/(exp(BUFFER[0]-BUFFER[1])+1);        
            Cavity[make_pair(i,j)][0] = 1 - x;
            Cavity[make_pair(i,j)][1] = x;           
/* 
            //normalize
            long double denom = 0;
	    long double md = 0; 
            for(int l =0;l<K;l++){
		md +=BUFFER[l];
	    }
	    md/=(long double) K;
            for(int l =0;l<K;l++){
		denom+= exp(BUFFER[l]-md);
	    }
            for(int l =0;l<K;l++){
                Cavity[make_pair(i,j)][l] = exp(BUFFER[l]-md) / denom;
	    }
*/
        }
    }
    return;
}
