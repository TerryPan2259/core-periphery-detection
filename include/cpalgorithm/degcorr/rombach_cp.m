classdef rombach_cp < cpabst 
	methods (Access = public)
		function [C,P,Q,Qs,score,param,dt] = detect( self, G, param )
			C= [];  P = [];Q = -Inf;score = [];
			dt = cputime;
		
			for it = 1:param.numRun
				if param.dcor  % degree correlation
					switch param.solver
						case 'label_switching'
							[Ct,Pt,Qt,scoret] = self.label_switching_dcor(G,param.alpha,param.beta);
						case 'annealing'
							[Ct,Pt,Qt,scoret] = self.simulated_annealing_dcor(G,param.alpha,param.beta);
						otherwise	
							disp('unknown solver for Rombach method')
					end
				else	
					switch param.solver
						case 'label_switching'
							[Ct,Pt,Qt,scoret] = self.label_switching(G,param.alpha,param.beta);
						case 'annealing'
							[Ct,Pt,Qt,scoret] = self.simulated_annealing(G,param.alpha,param.beta);
						otherwise	
							disp('unknown solver for Rombach method')
					end
				end
				
				if Qt > Q | it == 1
					C = Ct; P = Pt; Q = Qt; score = scoret; Qs = Qt;
				end
			end
			Qs = Q;
			dt = cputime-dt;
		end
		function [q,qs,score] = eval(self,G,x,varargin)
			q = x'*G.adjacency_matrix()*x;
			qs = q;
			score = x;
		end
		
		function param = initParam(self,param)
			if ~isfield(param,'disp');param.disp = false;end
			if ~isfield(param,'alpha');param.alpha = 1;end
			if ~isfield(param,'beta');param.beta = 0.5;end
			if ~isfield(param,'dcor');param.dcor = false;end
			if ~isfield(param,'solver');param.solver = 'label_switching';end
		
			if strcmp( param.name, 'rombach_dcor' )	
				param.dcor = true;
			end
		end
		
	end
	
	methods (Access = private)
		
		function [C,P,Q,score] = label_switching(self,G,alpha,beta)
			A = G.adjacency_matrix('binary');
			N = size(A,1);
			[r,c] = find(triu(A,1));
			score = rombach_label_switching( [ r,c ], N , length(r), alpha, beta );
			[val,ord] = sort(score,'descend');
			C = sparse(N,1,0);
			C( ord(1:ceil(N*(1-beta))) ) = 1;
			P = 1-C;
			Q = score'*A*score;
		end
		function [C,P,Q,score] = label_switching_dcor(self,G,alpha,beta)
			A = G.adjacency_matrix('binary');
			N = size(A,1);
			[r,c] = find(triu(A,1));
			score = rombach_label_switching_dcor( [ r,c ], N , length(r), alpha, beta );
			[val,ord] = sort(score,'descend');
			C = sparse(N,1,0);
			C( ord(1:ceil(N*(1-beta))) ) = 1;
			P = 1-C;
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
		
		function [C,P,Q,score] = simulated_annealing_dcor(self,G,alpha,beta)
			A = G.adjacency_matrix();
			deg = sum(A,2);
			M = sum(deg)/2;
			Mod = A -deg*deg'/(2*M)
			n = size(A,1);
			x = randsample(n,n);
			c = self.corevector(x,alpha,beta);	
				
			function x = gen(x)
				nids = randsample(n,2);
				x(nids) = x(nids(2:-1:1));	
			end	
			function q = fasteval(x)
				q = x'*Mod*x;
			end	
			clear aopt
			aopt.Generator = @gen;
			aopt.Verbosity = 0;
			[C,fval] = anneal(@(x)(-fasteval(x)),c,aopt);
			P = 1-C;
			score = C;
			Q = -fval;
		
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
