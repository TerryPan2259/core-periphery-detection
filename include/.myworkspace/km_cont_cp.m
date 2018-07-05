classdef km_cont_cp < cpabst
	
	methods ( Access = public )
		
		function [C, P, Q, Qs, score, param, cpu_time] = detect( self, G, varargin )
			% --------------------------------
			% Initialise
			% --------------------------------
			param = [];
			if nargin == 3
				param = varargin{1};
			end
			param = self.initParam( param );
			C = []; P = []; Qs = []; score = []; Q = - Inf;
			cpu_time = cputime;
	
			% --------------------------------
			% Optimise node labels using a label switching algorithm 
			% --------------------------------
			A = G.adjacency_matrix();
			for it = 1:param.numRun 
				ts = cputime;
				%[Ct, Pt, Qt, Qst,scoret] = self.one_cp_opt_algorithm( G, param );
				[Ct, Pt, Qt, Qst,scoret] = self.al_opt_algorithm( G, param ); 
				if Qt > Q
					C = Ct; P = Pt; Q = Qt; Qs = Qst;score = scoret; 
				end
				if param.disp; disp( sprintf( '%d/%d finished ( max Q = %f, %f seconds )', it, param.numRun, full( Q ), cputime - ts ) ); end
			end
			cpu_time = cputime - cpu_time;
		end
		
		function [Q, Qs, score] = eval( self, G, C, P, param )
			A = G.adjacency_matrix();
			deg = sum(A,2);
			M = sum(deg)/2;
	
			Qs = diag(C'* A * C)' - (deg'*C).^2 / (2*M);	
			Q = sum(Qs);
			score = sum(C,2);
		end
	
		function param = initParam( self, param )
			addpath /panfs/panasas01/emat/sk16722/gurobi/gurobi701/linux64/matlab/
			addpath(genpath('/panfs/panasas01/emat/sk16722/program/lib/YALMIP'));
			if isempty( param )	
				param.name = 'km';
			end	
			if ~isfield( param, 'numRun' ); param.numRun = 20; end
			if ~isfield( param, 'null_model' ); param.null_model = 'erdos-renyi'; end
			if ~isfield( param, 'disp' ); param.disp = false; end
		end
	end
	
	methods ( Access = private )
		% =========================	
		% KM algorithm for continuous version of core-periphery structure
		% null model = configuration model 
		% =========================	
		
		% Optimise c and x alternatively until no gain on Q is yield.	
		function [C, P, Q, Qs, score] = al_opt_algorithm( self, G, param )
	
			% read input netw
			A = G.adjacency_matrix();
			[r,c,v] = find(triu(A));
			L = [r,c,v];
			deg = G.degree();
			M = sum(deg)/2;
			N = size(A,1);
			
			% initialise
			x = deg./sqrt(deg);	
			x = x./sqrt(sum(x.^2));	
			C = sparse(1:N, 1:N, x);
			Q0 = self.eval( G, C, C, param ); 
	
			itnum = 0;
			while itnum <100 
	
				% optimise c	
				qmax = -Inf;
				Ctmp = km_config_cont_label_switching( L, N, size(L,1), 1:N, x);
				Ctmp = sparse(1:N,Ctmp, 1);
				
				% optimise x
				D = deg'*Ctmp;
				Ctmp = diag(deg) *Ctmp * diag(1./sqrt(D));
				Ctmp = Ctmp/sqrt(sum(Ctmp(:).^2));	
				
				% check convergence	
				Q1 = self.eval( G, Ctmp, Ctmp, param );
				if (Q1-Q0) < eps
					break;
				end
				C = Ctmp;
				Q0 = Q1;
				itnum = itnum + 1;
			end
			[Q,Qs,score] = self.eval( G, C, C, param ); 
			P = C;
			
		end	
		% Optimise c and x alternatively until no gain on Q is yield.	
		function [C, P, Q, Qs, score] = one_cp_opt_algorithm( self, G, param )
			A = G.adjacency_matrix();
			[r,c,v] = find(triu(A));
			L = [r,c,v];
			deg = G.degree();
			M = sum(deg)/2;
			N = size(A,1);
		
			% initialise
			itnum = 0;
			ns = N;
			As = A;
			ds = deg;
			C = [];
			R = sparse(1:N,1:N,1);
			while itnum <100
			
				% Find a core-periphery pair in networks
				xvar = sdpvar(ns,1);
				Constraints = [xvar'*xvar ==1, xvar>=0];
				J = (xvar'*As*xvar - sum(xvar.*ds).^2/(2*M));
				options = sdpsettings('verbose',0,'usex0',1);
				xs = rand(ns,1); xs = xs/sqrt(sum(xs.^2));
				assign(xvar,xs);
				sol = optimize(Constraints,-J,options);
				xs = value(xvar);
				s = (xs.^2 > 0.0001); 
				xs(~s) = 0;
				% formatting
				c = R * xs;
				C = [C,c];
				R = R(:,~s);
				As = A(~s,~s);
				ds = ds(~s);
				ns = sum(~s);
				if ns <1
					break;
				end
				itnum = itnum + 1;
			end
			[Q,Qs,score] = self.eval( G, C, C, param ); 
			P = C;
			
		end	
	end
	
end
