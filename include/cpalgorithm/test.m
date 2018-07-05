function test()
	
	addpath ../graph
	addpath ../randgraph
	addpath ../comalgorithm
	addpath ../sbm
	addpath(genpath('~/program/lib/YALMIP'));
	
	rg = randgraph();
	rgparam = rg.init('config');
	rgparam.type = 'soft';
	rgparam.numGraph = 1;	

	cpa = cpalgorithm();	
	km = cpa.init('km');
	km.null_model = 'config';
	km.numRun = 10;
	km.numRandGraph = 100;
	km.stat_test = 'qtest';
	km.pval = 0.05;
	%kmer.stat_test = 'oslom';
	%kmer.isSelfLoop = false;
	%rg = randgraph();
	%rgparam = rg.init('erdos-renyi');
	%rgparam.type = 'soft';	
	%kmer.rgparam = rgparam;
	
	G = graph( '~/program/data/poliblog/links.dat');	
	%G = graph( '~/program/data/karate/links.dat');	
	A = G.adjacency_matrix();
	A = A - diag(diag(A));
	A = sign(A + A');
	[C,P,Q,Qs] = cpa.detect(A, km);
	[slice,pvals,randGraphs,Crand,Prand] = cpa.statistical_test(C, P, A, km);
	pvals
	slice
end
