classdef fabio_cp 
	methods (Access = public)
		function [C,P,Q,Qs,score,model] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
			[C,P,Q,score] = self.greedy(graph(A));
			Qs = Q;
		end
		function q = eval(self,A,x)
			q = 1- (sum(x)-max(x))*2/(length(x)-2);
		end
	end
	methods (Access = private)
		function [C,P,Q,score] = greedy(self,G,model)
			n = G.sz();
			d = G.degree();
			m = G.edgenum();
			A = G.adjacency_matrix();
	
			
			
			x = zeros(n,1);
			[idx,val] = min(d);
			r = find( d==idx);
			idx = randsample(1:length(r),1);
			idx = r(idx);x(idx) = 1;
			ak = d(idx);bk = 0;
			alpha = zeros(n,1);
			alpha(1) = 0;
			
			for k = 2:n
				score = (2*ak * (x'*A)'  - bk * d)./(ak*(ak+d));	
				score(x>0) = Inf;	
				idx = min(score);	
				r = find( score ==idx );
				idx = randsample(1:length(r),1);
				idx = r(idx);
				x(idx) = 1;
				ak = ak + d(idx);bk = x'*A*x;
				alpha(idx) = bk/ak;
			end
			C = alpha;
			P = 1-alpha;
			Q = eval(self,A,alpha);	
			score = alpha;
		end
		function model = init_model(self,model)
			if ~isfield(model,'disp');model.disp = false;end
		end
	end
end
