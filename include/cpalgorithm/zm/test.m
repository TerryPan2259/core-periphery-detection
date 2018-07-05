function test()
	addpath ../../graph
	G = graph( '~/program/data/dolphin/links.dat');	
	%G = graph( '~/program/data/dolphin/links.dat');	
	A = G.adjacency_matrix('binary');

	[r,c,v] = find(triu(A));
	M = length(r);
	N = size(A,1)
		
	C = zm_bp( [r,c,v], N, M )	
	%bp(eList, C, N, 2, 50, 0.000001);
end
