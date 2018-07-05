// A C / C++ program for Dijkstra's single source shortest path algorithm.
// The program is for adjacency matrix representation of the graph
  
#include <algorithm>
#include <stdio.h>
#include <vector>
#include <limits.h>
#include <iostream>     // std::cout, std::end
#include <fstream>
  
// Number of vertices in the graph
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
  
// driver program to test above function
int main()
{
   /* Let us create the example graph discussed above */
   int graph[9][9] = {{0, 4, 0, 0, 0, 0, 0, 8, 0},
                      {4, 0, 8, 0, 0, 0, 0, 11, 0},
                      {0, 8, 0, 7, 0, 4, 0, 0, 2},
                      {0, 0, 7, 0, 9, 14, 0, 0, 0},
                      {0, 0, 0, 9, 0, 10, 0, 0, 0},
                      {0, 0, 4, 14, 10, 0, 2, 0, 0},
                      {0, 0, 0, 0, 0, 2, 0, 1, 6},
                      {8, 11, 0, 0, 0, 0, 1, 0, 7},
                      {0, 0, 2, 0, 0, 0, 6, 7, 0}
                     };

	
   std::vector<int> L;
   int N = 9;
   for(int i = 0;i < N;i++){
   	for(int j = 0;j < i;j++){
		if(graph[i][j]>0){
			int ind = MIN(i,j) + MAX(i,j) * N;
			L.push_back(ind);
		}
   	}	
   }	
   
     
    dijkstra(L, 0,N);
  
    return 0;
}
