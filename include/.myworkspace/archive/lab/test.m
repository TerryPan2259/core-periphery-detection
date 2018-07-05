function test()
	addpath ../../graph
	addpath ../../community
	addpath ../../random_graph
	addpath(genpath('/panfs/panasas01/emat/sk16722/program/lib/YALMIP'));

	k = 1;
	P =[
	1,1,
	1,0,
	];
	P = kron(eye(k),P);
	Nlist = kron(ones(1,k),[5,5]);
	P = P*1 + (1-P)*0;
	A = sbm(P,Nlist,1);	
	
	model.alpha = 0.8;model.beta = 0.5;	
	cpd = rombach_config();	
	[Ct,P,Q] = cpd.detect(A,model);
	Q
	cpd = rombach_cp();	
	[C,P,Q] = cpd.detect(A,model);
	full([C,Ct])
	Q
end
