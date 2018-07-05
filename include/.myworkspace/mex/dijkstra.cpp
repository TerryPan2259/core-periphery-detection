#include <algorithm>
#include <vector>
#include <set>
#include "mex.h"
#include "math.h"
#include <random>
#include <iostream>     // std::cout, std::end
#include <fstream>
#include <stdio.h>
#include <limits.h>
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

  
// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
int minDistance(std::vector<int> dist, std::vector<bool> sptSet)
{
   // Initialize min value
   int min = INT_MAX, min_index;
    int N = dist.size(); 
   for (int v = 0; v < N; v++)
     if (sptSet[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;
  
   return min_index;
}
  
// Funtion that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
std::vector<int> dijkstra(std::vector<int> graph, int src,int N)
{
     std::vector<int> dist;dist.assign(N,INT_MAX);     // The output array.  dist[i] will hold the shortest
                      // distance from src to i
     std::vector<bool>sptSet;sptSet.assign(N,false); // sptSet[i] will true if vertex i is included in shortest
                     // path tree or shortest distance from src to i is finalized

     // Distance of source vertex from itself is always 0
     dist[src] = 0;
  
     // Find shortest path for all vertices
     for (int count = 0; count < N-1; count++)
     {
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in first iteration.
       int u = minDistance(dist, sptSet);
  
       // Mark the picked vertex as processed
       sptSet[u] = true;
  
       // Update dist value of the adjacent vertices of the picked vertex.
       for (int v = 0; v < N; v++){
  
         // Update dist[v] only if is not in sptSet, there is an edge from 
         // u to v, and total weight of path from src to  v through u is 
         // smaller than current value of dist[v]
         
	int duv = 0;
	int ind = MIN(u,v) + MAX(u,v) * N;
	if(std::find(graph.begin(),graph.end(),ind)!=graph.end()){
		duv = 1; 
	}
		
         if (!sptSet[v] && duv && dist[u] != INT_MAX 
                                       && dist[u]+duv < dist[v])
            dist[v] = dist[u] + duv;
	}
     }
     return dist;
}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
   	double* edgeList =  mxGetPr(prhs[0]);
   	int N = (int) mxGetPr(prhs[1])[0];
   	int M = (int) mxGetPr(prhs[2])[0]; //length of edge list
	
		
	std::vector<int> L;
	for(int i = 0;i<M;i++){
		int rid = (edgeList[i]-1);
		int cid = (edgeList[i+M]-1);
		int ind = MIN(rid,cid) + MAX(rid,cid) * N;
		L.push_back(ind);
	}
    	
	plhs[0] = mxCreateDoubleMatrix( (mwSize)N, (mwSize)N, mxREAL); 
    	double* D = mxGetPr(plhs[0]);
	std::vector<int> dist;
	for(int i = 0;i< N;i++){
    		dist = dijkstra(L, i, N);
		for(int j = 0;j< N;j++){
			D[ i  +  j * N ] = dist[j];
		}
	}	
}
