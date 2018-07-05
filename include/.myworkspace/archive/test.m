function test()
	addpath ../graph
	addpath ../community
	addpath ../random_graph
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
		0,1,1,1,1,1,1,1;
		1,0,1,1,1,1,1,1;
		1,1,0,1,1,1,1,1;
		1,1,1,0,0,0,0,0;
		1,1,1,0,0,0,0,0;
		1,1,1,0,0,0,0,0;
		1,1,1,0,0,0,0,0;
		1,1,1,0,0,0,0,0
	];
	A = [
		0,1,1,1,1,1,1,1;
		1,0,1,1,1,1,1,1;
		1,1,0,0,0,0,0,0;
		1,1,0,0,0,0,0,0;
		1,1,0,0,0,0,0,0;
		1,1,0,0,0,0,0,0;
		1,1,0,0,0,0,0,0;
		1,1,0,0,0,0,0,0
	];

	k = 2;
	P =[
	1,1,
	1,0.1,
	];
	P = kron(eye(k),P);
	Nlist = kron(ones(1,k),[10,30]);
	P = P*1 + (1-P)*0;
	A = sbm(P,Nlist,1);	

	
	%T = readtable('~/program/data/poliblog/links.dat','Delimiter','\t');
	%G = graph(T);
	%A = G.adjacency_matrix();
	G = graph(A);
	M = G.modularity_matrix();
	cpd = cpdetection();
	com = community();
	clear model
	
	
	%model.name = 'borgatti';model.disp = true;
	%model.name = 'borgatti';model.solver = 'opt';model.disp = true;
	model.name = 'borgatti_multi';
	model.alpha = 1;
	%cmodel.name = 'modularity';
	tic
	model.beta  =0.5;
	[C,P,Q,score] = cpd.detect(A,model);
	toc
	full(C)
	pause
		
	denom=sqrt(sum(C.^2,1));
	C = C*diag(1./denom);
	trace(C'*M*C)	
	%full([(1:size(C,1))',C+P*2])	
end
