classdef borgatti_cp_diag 
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
			
			[C,P,Q,score] = self.Kernighan_Lin(A,model);
			
			Qs = Q;
		end
		function q = eval(self,A,x)
			n = size(A,1);
			cores = find(x);	
			peripheries = find(x==0);	
			L = [];
			if length(cores)>1
				t = nchoosek(cores,2);
				L =[L;t];
			end
			if length(peripheries)>1
				t = nchoosek(peripheries,2);
				L =[L;t];
			end
			
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
			for j = 1:n
				assigned =false(n,1);
				nc = sum(x);np = n-nc;
				cedges = x'*A*x/2;
				pedges = (1-x)'*A*(1-x)/2;
				Ec = cedges;	
				Ep = pedges;	
				pa = (Ec + Ep)/(nc*(nc-1)/2  + np*(np-1)/2);pb = (nc*(nc-1)/2)/(nc*(nc-1)/2  + np*(np-1)/2);
				Qold = (1-pa)*Ec - ((nc-1)*nc/2-Ec)*pa;
				Qold = Qold/( sqrt(pa*(1-pa)*pb*(1-pb)) );
				indeg = A*x;
				if isinf(Qold);Qold = 0;end;
				xnum = sum(x);	
	
				dQ = 0;dQmax = -Inf;
				for i = 1:ceil(n/2)
					% compute score
					score = sparse(n,1,0);
					nc = xnum +(1-2*x);np = n-nc;
					Ec = cedges + (1-2*x).*indeg;
					Ep = pedges + (1-2*(1-x)).*(deg-indeg);
					pa = (Ep+Ec)./(nc.*(nc-1)/2  + np.*(np-1)/2);
					pb = (nc.*(nc-1)/2)./(nc.*(nc-1)/2  + np.*(np-1)/2);
					score = (1-pa).*Ec - ((nc-1).*nc/2-Ec).*pa;
					score = score./( sqrt(pa.*(1-pa).*pb.*(1-pb)) );
					score(isinf(score)) = 0;
					score(isnan(score)) = 0;
				
					
					[qmax,idx] = max(score);
					assigned(idx) = true;
					cedges = cedges + (1-2*x(idx))*indeg(idx);
					pedges = pedges + (1-2*(1-x(idx)))*(deg(idx)-indeg(idx));
					indeg = indeg + (1-2*x(idx))*A(:,idx);
					xnum = xnum + (1-2*x(idx));
					x(idx) = 1-x(idx);
					dQ = dQ + qmax - Qold;
					Qold = qmax;
						
					if dQmax < dQ
						xbest = x;
						dQmax = dQ;
					end
				end
					
				if dQmax <=eps*10000;
					break;
				end
				x = xbest;
				C = xbest;
			end
			P = 1-C;Q = self.eval(A,C);score = C;	
		end

		
		
		function model = init_model(self,model)
			if ~isfield(model,'solver');model.solver = 'kl';end
			if ~isfield(model,'loopnum');model.loopnum = 10000;end
		end
		
	end
end
