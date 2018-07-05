function test()
	addpath ../../graph
	addpath ../../sbm
	addpath ../../randgraph
	addpath(genpath('~/program/lib/YALMIP'));

	A = B;
	%C = [A,zeros(size(A));zeros(size(B)),B];
	A = kron(eye(2),B);
	N = size(A,1);

	
	%G = graph('~/program/data/karate/links.dat');
	%A = G.adjacency_matrix();
	
	cpa = maxcut_cp();
	param.name = 'maxcut';
	param = cpa.initParam(param);
	param.numRun = 1;
	param.optalgorithm='kl';
	param.temperature = 10;
		
		
	G = graph(A);	
	T = G.transition_matrix();	
	s = G.stationary_distribution();
	W  = T * diag(s);
	[C,P,Q] = cpa.detect( graph(A), param );
	full(C+P*2)
end
