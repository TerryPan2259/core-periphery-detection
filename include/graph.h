#include <iostream>
#include <vector>

using namespace std;
//typedef pair<int, double> Edge; // Edge 

class Neighbour{
public:
    int _node;
    double _w; 
    // constracter
    Neighbour();
    Neighbour(int node, double w);
	
    // Getters
    int get_node() const;
    double get_w() const;
};

Neighbour::Neighbour(int node, double w){
	_node = node;
	_w = w;
};

// Getter -----------
int Neighbour::get_node() const{
	return _node;
}

double Neighbour::get_w() const{
	return _w;
}

    
class Graph{
public:
    vector<vector<Neighbour>> _neighbours;
    
    // constracter
    Graph();
    
    // constracter
    Graph(int num_nodes);
	
    // Getters
    int get_num_nodes() const;
    int get_num_edges() const;
    int degree(int nid) const;
    double wdegree(int nid) const;
    void get_weight(int nid, int j, int& nei, double& w) const;
    
    void compress();

    Neighbour get_kth_neighbour(int nid, int k) const;

    void print() const;
    vector<vector<double>> to_matrix() const;
    
    void addEdge(int u, int v, double w);
    
};

Graph::Graph(){
	vector<vector<Neighbour>> tmp(0, vector<Neighbour>(0));	
	_neighbours = tmp;
};
    
Graph::Graph(int num_nodes){
	vector<vector<Neighbour>> tmp(num_nodes, vector<Neighbour>(0));	
	_neighbours = tmp;
};


// Getter -----------
int Graph::get_num_nodes() const{
	return _neighbours.size();
}

int Graph::get_num_edges() const{
	int N = get_num_nodes();
	int M = 0;
    	for (int i = 0; i < N; i++) {
		int sz = degree(i);
    		for (int j = 0; j < sz; j++) {
			M+= _neighbours[i][j].get_w();
		}
	}
	return M/2;
}

// get the id and weight of the jth neighbour of node nid
void Graph::get_weight(int nid, int j, int& nei, double& w) const{
	nei = _neighbours[nid][j].get_node();
	w = _neighbours[nid][j].get_w();
}

// get weighted degree
double Graph::wdegree(int nid) const{
	int sz = _neighbours[nid].size();
	double deg = 0;
    	for (int j = 0; j < sz; j++) {
        	deg+= _neighbours[nid][j].get_w();
	}
	return deg;
}

// get degree
int Graph::degree(int nid) const{
	return _neighbours[nid].size();
}

Neighbour Graph::get_kth_neighbour(int nid, int k) const{
	return _neighbours[nid][k];
}

// add 
void Graph::addEdge(int u, int v, double w){
	
	int sz = _neighbours.size();
	while( (sz <=u) ){
		vector<Neighbour> tmp(0);
		_neighbours.push_back(tmp);
		sz++;
	}

	Neighbour ed1(v, w);
	_neighbours[u].push_back(ed1);

	/*	
	if(u==v){
		Neighbour ed1(v, w);
		_neighbours[u].push_back(ed1);
	}else{
		Neighbour ed1(v, w);
		_neighbours[u].push_back(ed1);

		Neighbour ed2(u, w);
		_neighbours[v].push_back(ed2);
	}*/
}

// add 
void Graph::print() const{
	int N = get_num_nodes();
	for(int i =0; i < N;i++){
		int sz = _neighbours[i].size();
		for(int j =0; j < sz;j++){
			cout<<i<<" "<<_neighbours[i][j].get_node()<<" "<<_neighbours[i][j].get_w()<<endl;
		}
	
	}
}

// add 
vector<vector<double>> Graph::to_matrix() const{
	int N = get_num_nodes();
	vector<vector<double>> M(N, vector<double>(N, 0));
		
	for(int i =0; i < N;i++){
		int sz = _neighbours[i].size();
		for(int j =0; j < sz;j++){	
			M[i][_neighbours[i][j].get_node()] = _neighbours[i][j].get_w();
		}
	}
	return M;
}

// merge multiple-edges with a single weighted edge 
void Graph::compress(){
	
	int N = get_num_nodes();
	
	// copy	
        vector<vector<Neighbour>> prev_neighbours = _neighbours;

	// initialise _neighbours	
	for(int i = 0; i < N; i++){
		_neighbours[i].clear();
	}
	_neighbours.clear();
	vector<vector<Neighbour>> tmp(N, vector<Neighbour>(0));	
	_neighbours = tmp;

	for(int i = 0; i < N; i ++){
		int sz = prev_neighbours[i].size();
		map<int, double> myMap;
		for(int j =0; j < sz;j++){	
			int nei = prev_neighbours[i][j].get_node();
			double w = prev_neighbours[i][j].get_w();
			if ( !myMap.insert( make_pair( nei, w ) ).second ) {
				myMap[nei]+=w;
			}
		}
		for (const auto & p : myMap) {
			addEdge(i, p.first, p.second);	
		}
	}
}
