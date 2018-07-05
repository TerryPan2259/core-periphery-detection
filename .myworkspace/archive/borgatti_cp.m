classdef borgatti_cp 
	methods (Access = public)
		function [C,P,Q,Qs,score,model] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
			C= [];  P = [];Q = [];score = [];
			switch model.solver
				case 'mcmc'
					[C,P,Q,score] = self.mcmc(graph(A),model);
				case 'kl'
					[C,P,Q,score] = self.Kernighan_Lin(A,model);
				case 'eigen'
					[C,P,Q,score] = self.eigen(graph(A));
				case 'anealing'
					[C,P,Q,score] = self.simulated_anealing(A);
				case 'ga'
					[C,P,Q,score] = self.genetic_algorithm(A);
				case 'opt'
					[C,P,Q,score] = self.opt(graph(A));
				otherwise	
					disp('unknown solver for borgatt & Everett model')
			end
			Qs = Q;
		end
		function q = eval(self,A,x)
			n = size(A,1);
			L = nchoosek(1:n,2);
			eids = sub2ind([n,n],L(:,1),L(:,2));
			a = A(eids);	
			xx = x(L(:,1)) + x(L(:,2)) - ( x(L(:,1)) .* x(L(:,2)) );	
			q = corr(a,xx);	
			
			if isnan(q)
				q = 0;
			end	
		end
	end
	methods (Access = private)
		function [C,P,Q,score] = Kernighan_Lin(self,A,model)
			n = size(A,1);
			deg = sum(A,2);
			M = n*(n-1)/2;
			muy = sum(deg)/2;
			
			Q = 0;
			% init x
			if n<=2
				C = ones(n,1);P=1-C;Q=1;score = zeros(n,1);
				return;
			end

			while(true)
				x = round(rand(n,1));
				if ~all(x) & ~all(~x);break;end
			end
		
			
			C = x;
			% update x
			nc = sum(x);
			m = nc *(nc-1)/2 + nc * (n-nc);
			indeg = A*x;
			dot = x'*deg - x'*indeg/2; 
			for j = 1:n	
				assigned =false(n,1);
				Qold = (dot - m * muy/M)/sqrt(m *(M-m));
				dQ = 0;dQmax = -Inf;
				for i = 1:ceil(n/2)
					dm = (nc-n).*x + (n-nc-1).*(1-x);
					ddot = (indeg-deg); ddot(x==0) = -ddot(x==0);
					ev = (dot+ddot - (m+dm) * muy/M) ./ sqrt( (m+dm).*(M-(m+dm)));
					ev(assigned)=-Inf;ev(isinf(ev)) = -Inf;
					
					[qmax,idx] = max(ev);
					assigned(idx) = true;
					
					nc = nc - (2*x(idx)-1);
					m = m + dm(idx);				
					dot = dot + ddot(idx);
					indeg = indeg + A(:,idx) * (-2*x(idx)+1);	
	
					x(idx) = 1-x(idx);
					dQ = dQ + qmax - Qold;
					Qold = qmax;
					
					if dQmax < dQ
						xbest = x;
						dQmax = dQ;
						ncbest = nc;mbest = m;dotbest = dot;indegbest = indeg;
					end
				end
					
				if dQmax <=eps*10000;
					break;
				end
				nc = ncbest;m = mbest;dot = dotbest;indeg = indegbest;
				x = xbest;
				C = xbest;
			end
			P = 1-C;Q = self.eval(A,C);score = C;	
		end

		function [C,P,Q,score] = genetic_algorithm(self,A,model)
			n = size(A,1);
			L = nchoosek(1:n,2);
			eids = sub2ind([n,n],L(:,1),L(:,2));
			a = A(eids);	
			function q = evalcp(x)
				xx = any([x(L(:,1));x(L(:,2))]',2);	
				q = -corr(a,xx);	
			end
			
			clear problem;	
			problem.fitnessfcn = @evalcp;
			problem.nvars = n;
			problem.lb = zeros(n,1);
			problem.ub = ones(n,1);
			problem.intcon = 1:n;
			problem.solver ='ga';
			problem.options = optimoptions('ga','UseParallel', true, 'UseVectorized', false);
			C = ga(problem);
			P = 1-C;
			Q = evalcp(C);
			C = C';P=P';
			score = sum(A,1)';	
		end
		function [C,P,Q,score] = simulated_anealing(self,A,model)
			n = size(A,1);
			L = nchoosek(1:n,2);
			eids = sub2ind([n,n],L(:,1),L(:,2));
			a = A(eids);	
			function q = evalcp(x)
				xx = any([x(L(:,1)),x(L(:,2))],2);	
				q = -corr(a,xx);	
			end
			
			xinit = self.eigen(graph(A));	
			C = simulannealbnd(@evalcp,xinit);
			P = 1-C;
			Q = evalcp(C);
			score = sum(A,1)';	
		end

		function [C,P,Q,score] = mcmc(self,G,model)
			n = G.sz();	
			A = G.adjacency_matrix();
			
			C = rand(n,1)<0.5;
			p = 0;
			
			itnum = 0;
			L = nchoosek(1:n,2);
			eids = sub2ind([n,n],L(:,1),L(:,2));
			Cmax = C;pmax = p;
			B = ones(n);a = A(eids);	
			while(itnum < model.loopnum)
				
				nid = randsample(n,1);
				
				C(nid) = 1-C(nid);
				X = B;
				X(C==0, C==0) = 0;
				
				ptest = corr(a,X(eids));	
					
				if rand() < (ptest+1)/(p+1) % add one to prevent p to be minus
					p = ptest;
					if pmax < p
						Cmax = C;
						pmax = p;
					end
				else
					C(nid) = 1-C(nid);
				end	
				itnum = itnum + 1;
			end
			C = Cmax;
			P = 1-Cmax;
			Q = pmax;
			score = sum(A,1)';	
		end
		
		function [C,P,Q,score] = opt(self,G)
			deg = G.degree();	
			
			[val,ord] = sort(deg,'descend');
			dsum = 0;
			n = G.sz();p = G.density();A = G.adjacency_matrix();
			q = zeros(length(ord),1);
			C = sparse(length(ord),1,0);
			for nc = 1:length(ord)-2
				dsum = dsum + val(nc);
				C(ord(nc)) = 1;
				given = dsum - C'*A*C/2;
				null = nc*(n-nc) + nc*(nc-1)/2;
				
				vy = null / (n*(n-1)/2);
				q(nc) = (given - p*null)/ sqrt(vy*(1-vy));
			end
			[~,idx] = max(q);
			C = sparse(ord(1:idx),1,1,n,1);
			P = 1-C;
			Q = self.eval(G.adjacency_matrix(),C); 
			score = C;
		end		
	
		function [C,P,Q,score] = eigen(self,G)
			cmp = G.connected_components();
			N = G.N;
			% compute minimum eigen vector of modularity matrix for each components 
			Q = [];
			C = [];P=[];	
			for c = cmp 
				At= G.A(c>0,c>0);
				v = maximum_eigen(At);
				v = abs(v);
				
				[vals,idx] = sort(v,'descend');
				
				xbest = sparse(length(v),1,0);
				bestscore = 0;	
				for i = 1:length(v)
					xcand = sparse(length(v),1,0);xcand(v>=vals(i))=1;
					template = xcand*xcand'+ xcand*(1-xcand)' + (1-xcand)*xcand';
					template = template - diag(diag(template));
					score = corr(At(:),template(:));
					if score >= bestscore
						bestscore = score;
						xbest = xcand;
					end	
				end	
				tmp = sparse(N,1,0);
				tmp(c>0) = xbest;
					
				C = [C,tmp];
				
				tmp = sparse(N,1,0);
				tmp(c>0) = 1-xbest;
				P = [P,tmp];
				Q = [Q, bestscore];
			end
			C = sum(C,2);
			P = 1-C;
			
			if length(C)==0
				Q = 0;
				C = sparse(G.sz(),1,0);
				P = 1-C;score = C;
				return;
			end
			Q = max(Q);
			score = sum(G.A,2);
		end
		
		
		function model = init_model(self,model)
			if ~isfield(model,'solver');model.solver = 'kl';end
			if ~isfield(model,'loopnum');model.loopnum = 10000;end
		end
		
	end
end
