classdef rombach_cp < cpabst 
	methods (Access = public)
		function [C,P,Q,Qs,score,param,dt] = detect( self, G, param )
			C= [];  P = [];Q = -Inf;score = [];
			dt = cputime;
		
			for it = 1:param.numRun
				if strcmp( param.null_model, 'config')  % degree correlation
					switch param.solver
						case 'label_switching'
							[Ct,Pt,Qt,scoret] = self.label_switching_dcor(G,param);
							Qst = Qt;
						case 'annealing'
							[Ct,Pt,Qt, Qst, scoret] = self.simulated_annealing_dcor(G,param);
						otherwise	
							disp('unknown solver for Rombach method')
					end
				elseif strcmp( param.null_model, 'erdos-renyi')	
					switch param.solver
						case 'label_switching'
							[Ct,Pt,Qt,scoret] = self.label_switching(G,param.alpha,param.beta);
							Qst = Qt;
						case 'annealing'
							[Ct,Pt,Qt,scoret] = self.simulated_annealing(G,param.alpha,param.beta);
							Qst = Qt;
						otherwise	
							disp('unknown solver for Rombach method')
					end
				end
				
				if Qt > Q | it == 1
					C = Ct; P = Pt; Q = Qt; score = scoret; Qs = Qst;
				end
			end
			Qs = Q;
			dt = cputime-dt;
		end
		function [Q,Qs,score] = eval(self,G,C,P,param)
			A = G.adjacency_matrix('binary');
			if strcmp( param.null_model, 'erdos-renyi')
				Q = C'*A*C;
				Qs = Q;
			elseif strcmp( param.null_model, 'config')
				deg = sum(A,2);
				M = sum(deg)/2;
				Qs = diag(C' *A * C) - (C'*deg).^2/(2*M);
				Q = sum(Qs);
			end
			score = sum(C,2);
		end
		
		function param = initParam(self,param)
			if ~isfield(param,'disp');param.disp = false;end
			if ~isfield(param,'alpha');param.alpha = 0.5;end
			if ~isfield(param,'beta');param.beta = 0.8;end
			if ~isfield(param,'null_model');param.null_model = 'erdos-renyi';end
			if ~isfield(param,'numRun');param.numRun = 10;end
			if ~isfield(param,'solver');param.solver = 'annealing';end
		end
		
	end
	
	methods (Access = private)
		
		function [C,P,Q,score] = label_switching(self,G,alpha,beta)
			A = G.adjacency_matrix('binary');
			N = size(A,1);
			[r,c] = find(triu(A,1));
			score = rombach_label_switching( [ r,c ], N , length(r), alpha, beta );
			C = score; P = 1-score;
			Q = score'*A*score;
		end
		function [C,P,Q,score] = label_switching_dcor(self,G,alpha,beta)
			A = G.adjacency_matrix('binary');
			N = size(A,1);
			[r,c] = find(triu(A,1));
			score = rombach_label_switching_dcor( [ r,c ], N , length(r), alpha, beta );
			C = score; P = 1-score;
			deg = G.degree();
			M = G.numEdge();	
			Q = score'*A*score-sum(score.*deg).^2/(2*M);
		end
		
		
		function [C,P,Q,score] = simulated_annealing(self,G,alpha,beta)
			A = G.adjacency_matrix();
			n = size(A,1);
			x = randsample(n,n);
			c = self.corevector(x,alpha,beta);	
				
			function x = gen(x)
				nids = randsample(n,2);
				x(nids) = x(nids(2:-1:1));	
			end	
			function q = fasteval(x)
				q = x'*A*x;
			end	
			clear aopt
			aopt.Generator = @gen;
			aopt.Verbosity = 0;
			[C,fval] = anneal(@(x)(-fasteval(x)),c,aopt);
			P = 1-C;
			score = C;
			Q = -fval;
		
		end
		
		function [C,P,Q,Qs,score] = simulated_annealing_dcor(self,G,param)
			alpha = param.alpha;	
			beta = param.beta;
			A = G.adjacency_matrix();
			[r,c,v] = find(triu(A));
			L = [r,c,v];
			deg = sum(A,2);
			M = sum(deg)/2;
			N = size(A,1);
			x = randsample(N,N);
			coreness = self.corevector(x,alpha,beta);	
				
			function x = gen(x)
				nids = randsample(N,2);
				x(nids) = x(nids(2:-1:1));	
			end	
			function q = fasteval(x,C)
				q = trace( C'*diag(x) * A * diag(x)*C ) - sum(C'*(deg.*x).^2)/(2*M);
			end	
			clear aopt
			aopt.Generator = @gen;
			aopt.Verbosity = 0;
		
			C = sparse(1:N,1:N,coreness);	
			itnum = 0;
			while itnum < 100
				Q0 = self.eval( G, C, C, param );
	
				ctmp = km_config_cont_label_switching( L, N, size(L,1), 1:N, coreness);
				Ctmp = sparse(1:N,ctmp,1);
				[coreness_tmp,fval] = anneal(@(x)(-fasteval(x,Ctmp)),coreness,aopt);
			
				% check convergence
				Ctmp = diag(coreness_tmp) * Ctmp;	
				Q1 = self.eval( G, Ctmp, Ctmp, param );
				[Q0,Q1]
				if (Q1-Q0) < eps
					break;
				end
				C = Ctmp;
				coreness = sum(C,2);
				itnum = itnum + 1;
					
			end
			P = C;	
			[Q,Qs,score] = self.eval( G, C, C, param );
		end
		
		
		function c = corevector(self,x,alpha,beta)
			n = length(x);
			bn = floor(beta*n);
			cx = x<=bn;
			px = ~cx;
			c = (1-alpha)/(2*bn)*x.*cx + ((x.*px-bn)*(1-alpha)/(2*(n-bn))  + (1+alpha)/2).*px;
		end
					
	end
end
