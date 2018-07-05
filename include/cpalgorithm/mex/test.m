	
	addpath ../../sbm
	addpath ../../graph
	addpath ../
		
	G = graph('~/program/data/poliblog/links.dat');	
	A = G.adjacency_matrix('binary');

	A = G.adjacency_matrix();	
	[r,c,v] = find(triu(A,1));
	dt = cputime;
	C = bealgorithm( [r,c], size(A,1), length(r) );
	boyd_stest( [r,c], size(A,1), length(r), C, 0.01, 100 )
	dt = cputime-dt;
	disp(dt)
