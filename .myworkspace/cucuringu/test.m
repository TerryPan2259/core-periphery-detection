	
	addpath ../../graph
	T = readtable('~/program/data/karate/links.dat','Delimiter','\t');
	G = graph(T);	
	A = G.adjacency_matrix();
	
	A = ones(100);
	A(51:end,51:end) = 0;
	%A = ones(10);
	%A(6:end,6:end) = 0;
	
	A = A - diag(diag(A));	
	G = graph(A);	
	ccp = cucuringu_cp();
	param.name ='low_rank_core';
	param = ccp.initParam(param);
	tic	
	C = ccp.detect(G,param);	
	toc
	param.name ='lap_core';
	param = ccp.initParam(param);
	tic	
	C = ccp.detect(G,param);	
	toc
	param.name ='lapsgn_core';
	param = ccp.initParam(param);
	tic	
	C = ccp.detect(G,param);	
	toc
	%C = dijkstra([r,c],N,length(r));
	%C = minres_discrete([r,c],N,length(r));
