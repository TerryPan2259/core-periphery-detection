Included are: 
Belief Propagation algorithm implemented in C++. 
Python script for converting 'gml' format to 'edgelist' format to make it readable for C++.
Things to be fixed manuelly in the C++ program:
N: number of nodes of the network, copy from the network info of python script.
M: number of edges of the network, copy from the newtork info of python script.
inputfile: edgelist file for the network.
outputfile: output group assignment list for the core/periphery division.

C++11 is required to compile. In ubuntu or linux compile with:
g++ -std=c++0x bp.cpp


