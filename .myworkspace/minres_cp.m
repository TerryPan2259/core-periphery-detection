classdef minres_cp < cpabst 
	methods (Access = public)
		function [C,P,Q,Qs,score,param,dt] = detect(self, G, param)
			C= [];  P = [];Q = Inf;score = [];
			ts = cputime;
			
			for it = 1:param.numRun
				switch param.name
					case 'minres_cont'
						[Ct,Pt,Qt,Qst,scoret] = self.grad_descent(G,param); 
					case 'minres'
						[Ct,Pt,Qt,Qst,scoret] = self.Kernighan_Lin(G,param); 
					otherwise	
				end
				if Qt < Q % note that the smaller is better
					C = Ct;P = Pt;Q = Qt;score = scoret;Qs = Qst;
				end
			end
			dt = cputime-ts;
		end

		function [Q,Qs,score] = eval(self,G,x,varargin)
			name = 'minres';
			if nargin == 5
				param = varargin{2};
				if isfield(param,'name');name = param.name;end
			end
			if strcmp(name, 'minres_cont')	
				A = G.adjacency_matrix();
				M = G.numEdge();
				Q = 2*M-2*x'*A*x + (x'*x)^2 - sum(x.^4);
				Qs = Q;score = x;
			else 
				A = G.adjacency_matrix('binary');
				mcc= x'*A*x;
				mpp= (1-x)'*A*(1-x);
				nc = sum(x);
				np = size(A,1)-nc;
			
				Q = (nc*(nc-1) -mcc) + mpp;
				Qs = Q; score = sum(A,2);
			end
			Q = -Q;Qs = -Qs;score = -score; 
		end
		
		function param = initParam(self,param)
			if ~isfield(param,'name');param.solver = 'minres';end
			if ~isfield(param,'numRun');param.numRun = 10;end
			if ~isfield(param,'disp');param.disp = false;end
		end
		
	end
	methods (Access = private)
		function [C,P,Q,Qs,score] = Kernighan_Lin(self,G,param)
			A = G.adjacency_matrix('binary');
			deg = sum(A,2);
			N = size(A,1);	
			[val,nids] = sort(deg,'descend');
			Z = sum(deg)/2;Zbest = Inf;kbest = 0;
			for k = 1:N-1
				Z = Z + k - 1 - val(k);
				if Z < Zbest
					kbest = k;
					Zbest = Z;	
				end	
			end
			C = sparse(nids(1:kbest),1,1,N,1);
			P = 1-C;
			[Q,Qs,score] = self.eval( G, C, P, param); 
		end
		
		function [C,P,Q,Qs,score]=minressvd(self,G,param)
			A = G.adjacency_matrix('binary');
			C = []; P =[];Q = []; Qs = [];score = [];
			n = size(A,1);
			u = sdpvar(n,1);	
			Constraints = [u>=0];
				
			Objective = (sum(sum((A- u*u').^2)) - sum((diag(A)-u.*u).^2)); 
			options = sdpsettings('verbose',0);
			while(true)
				sol = optimize(Constraints,Objective,options);
				% Analyze error flags
				if sol.problem == 0
					% Extract and display value
					u = value(u);
					break
				else
					display('Hmm, something went wrong!');
					sol.info
					yalmiperror(sol.problem)
				end
			end
			d = sqrt(u'*u)*sqrt(u'*u);
			C = u;
			P = 1-C;
			Q = (sum(sum((A- u*u').^2)) - sum((diag(A)-u.*u).^2)); 
			Qs = Q;
			score  = u;	
			Q = -Q;Qs = -Qs;score = -score; 
		end
%		function [C,P,Q,Qs,score] = Kernighan_Lin(self,G,param)
%			A = G.adjacency_matrix('binary');
%			n = size(A,1);
%			deg = sum(A,2);
%			M = n*(n-1)/2;
%			muy = sum(deg)/2;
%			
%			Q = 0;
%			% init x
%			if n<=2
%				C = ones(n,1);P=1-C;Q=1;score = zeros(n,1);
%				return;
%			end
%
%			while(true)
%				x = round(rand(n,1));
%				if ~all(x) & ~all(~x);break;end
%			end
%		
%			
%			C = x;
%			% update x
%			nc = sum(x);
%			m = x'*A*x;
%			tocore = A*x;
%			for j = 1:n	
%				assigned =false(n,1);
%				Qold = -self.eval(G,x);
%				dQ = 0;dQmax = -Inf;
%				for i = 1:ceil(n/2)
%					
%					toperi = deg-tocore;
%				
%					ev = -((2*x-1).*(tocore + toperi - nc) + x);	
%					
%					ev(assigned)=-Inf;ev(isinf(ev)) = -Inf;
%					
%					[qmax,idx] = max(ev);
%					assigned(idx) = true;
%					
%					nc = nc - (2*x(idx)-1);
%					tocore = tocore + A(:,idx) * (-2*x(idx)+1);	
%	
%					x(idx) = 1-x(idx);
%					dQ = dQ + qmax - Qold;
%					Qold = qmax;
%					
%					if dQmax < dQ
%						xbest = x;
%						dQmax = dQ;
%						ncbest = nc;tocorebest = tocore;
%					end
%				end
%				if dQmax <=eps*10000;
%					break;
%				end
%				nc = ncbest;tocore = tocorebest;
%				x = xbest;
%				C = xbest;
%			end
%			P = 1-C;[Q,Qs,score] = self.eval(G,C,P);	
%		end
		
	end
end
