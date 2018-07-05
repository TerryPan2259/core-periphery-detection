classdef minres_cp < cpabst 
	methods (Access = public)
		function [C,P,Q,Qs,score,param,dt] = detect(self, G, param)
			C= [];  P = [];Q = Inf;score = [];
			ts = cputime;
			
			if ~isfield(param,'dcor'); param.dcor = false; end
				
			if param.dcor 
				for it = 1:param.numRun
					switch param.name
						case 'minres_cont_dcor'
							[Ct,Pt,Qt,Qst,scoret] = self.grad_descent_dcor(G,param); 
						otherwise	
					end
					if Qt < Q % note that the smaller is better
						C = Ct;P = Pt;Q = Qt;score = scoret;Qs = Qst;
					end
				end
			else
				for it = 1:param.numRun
					switch param.name
						case 'minres_cont'
							[Ct,Pt,Qt,Qst,scoret] = self.grad_descent(G,param); 
						case 'minres'
							[Ct,Pt,Qt,Qst,scoret] = self.Kernighan_Lin(G,param); 
						case 'minres_discrete'
							[Ct,Pt,Qt,Qst,scoret] = self.Kernighan_Lin(G,param);
						otherwise	
					end
					if Qt < Q % note that the smaller is better
						C = Ct;P = Pt;Q = Qt;score = scoret;Qs = Qst;
					end
				end
			end
			dt = cputime-ts;
		end

		function [q,qs,score] = eval(self,G,x,varargin)
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
			elseif strcmp(name, 'minres_cont_dcor')	
				A = G.adjacency_matrix();
				M = G.numEdge();
				deg = G.degree();	
				Q = (-x'*A*x + ((deg'*x)^2)/(2*M));
				Qs = Q;score = x;
			else 
				A = G.adjacency_matrix('binary');
				mcc= x'*A*x;
				mpp= (1-x)'*A*(1-x);
				nc = sum(x);
				np = size(A,1)-nc;
			
				q = (nc*(nc-1) -mcc) + mpp;
				qs = q; score = sum(A,2);
			end 
		end
		
		function param = initParam(self,param)
			if ~isfield(param,'name');param.solver = 'minres';end
			if ~isfield(param,'numRun');param.numRun = 10;end
			if ~isfield(param,'dcor');param.dcor = false;end
			
			% for grad descent 
			if ~isfield(param,'sigma');param.sigma = 0.1;end
			if ~isfield(param,'beta');param.beta = 0.1;end
			if ~isfield(param,'mu');param.mu = 1e-4;end
			if strcmp( param.name, 'minres_cont_dcor' )	
				param.dcor = true;
			end
			
		end
		
	end
	methods (Access = private)
		
		function [C,P,Q,Qs,score] =  grad_descent(self, G, param)
			
			A = G.adjacency_matrix('binary');
			N = size(A,1);	
			deg = G.degree(); M = G.numEdge();
			
			x = rand(N,1);
			%Qx = trace( (A-x*x')' * (A-x*x') ); 
			Qx = 2*M-2*x'*A*x + (x'*x)^2 - sum(x.^4);
			while(true)
				alpha = 1;	
				while(true)
					%S = x'*x * eye(N) - diag(x.^2);
					%df = param.sigma * (-A* x + 2*S*x);
					df = param.sigma * (-A* x + (x'*x)*x - (x.^2).*x);
					
					xnew = x - alpha * df;
					xnew(xnew<0) = 0;
					%xnew = xnew/sqrt(xnew'*xnew);
					%Qxnew = trace( (A-xnew*xnew')' * (A-xnew*xnew') ); 
					Qxnew = 2*M-2*xnew'*A*xnew + (xnew'*xnew)^2 - sum(xnew.^4);
					if Qxnew  - Qx  <= param.sigma * df'*(xnew-x)
						break;
					end
					alpha = alpha * param.beta; 
				end	
				if abs(Qxnew  - Qx)/N < param.mu; 
					x = xnew;
					Qx = Qxnew;
					break;
				end
				x = xnew;
				Qx = Qxnew;
			end
			C = x;
			P = 1-C;	
			Q = Qx;
			Qs = Q;
			score = x;	
		end
		
		function [C, P, Q, Qs, score] = grad_descent_dcor( self, G, param )
			
			A = G.adjacency_matrix('binary');
			N = size(A,1);	
			deg = G.degree(); M = G.numEdge();
			
			x = rand(N,1);
			%Zx = x/sqrt(x'*x);
			Qx = (-x'*A*x + ((deg'*x)^2)/(2*M));
			while(true)
				alpha = 1;	
				while(true)
					df = param.sigma * (-A* x + (deg)*(deg'*x)/(2*M));
					%df = sigma * (Mod * x);
					xnew = x - alpha * df;
					xnew(xnew<0) = 0;
					%xnew = xnew/sqrt(xnew'*xnew);
					Qxnew = (-xnew'*A*xnew + ((deg'*xnew)^2)/(2*M));
					if Qxnew  - Qx  <= param.sigma * df'*(xnew-x)
						break;
					end
					alpha = alpha * param.beta; 
				end
				if abs(Qxnew  - Qx)/N < param.mu; 
					x = xnew;
					Qx = Qxnew;
					break;
				end
				x = xnew;
				Qx = Qxnew;
			end
			C = x;
			P = 1-C;	
			Q = -Qx;
			Qs = Q;
			score = x;	
			
		end
		
		
		
		function [C,P,Q,Qs,score] = Kernighan_Lin(self,G,param)
			A = G.adjacency_matrix('binary');
			[r,c] = find(triu(A,1));
			C = minres_discrete( [r,c], size(A,1), length(r) );
			P = 1-C;[Q,Qs,score] = self.eval(G,C,P);	
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
		end
		

	end
end
