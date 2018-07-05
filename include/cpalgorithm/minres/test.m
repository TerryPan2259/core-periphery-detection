function test()
	addpath ../../graph
	addpath(genpath('~/program/lib/YALMIP'));
	addpath(genpath('/panfs/panasas01/emat/sk16722/program/lib/YALMIP'));

	A = [
		0,1,1,1,1,0,0,0;
		1,0,1,1,0,1,0,0;
		1,1,0,1,0,0,1,0;
		1,1,1,0,0,0,0,1;
		1,0,0,0,0,0,0,0;
		0,1,0,0,0,0,0,0;
		0,0,1,0,0,0,0,0;
		0,0,0,1,0,0,0,0
	];
	A = [
		0,1,1,1,1,1;
		1,0,1,1,1,1;
		1,1,0,1,1,1;
		1,1,1,0,0,0;
		1,1,1,0,0,0;
		1,1,1,0,0,0;
	];
	A = kron(eye(2),A);
	%T = readtable('~/program/data/karate/links.dat','Delimiter','\t');
	%G = graph(T);	
	%A = G.adjacency_matrix();
	
	cpa = minres_cp();
	param.name = 'minres_cont';
	param = cpa.initParam(param);
	param.numRun = 1;
	
	[C,P,Q,Qs,score,param,cpu_time] = cpa.detect(graph(A),param);	
	C
	score
	Q
end
