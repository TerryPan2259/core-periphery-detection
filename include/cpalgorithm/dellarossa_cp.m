classdef dellarossa_cp < cpabst 
	methods (Access = public)
		function [C,P,Q,Qs,score,param,dt] = detect(self, G, param)
			dt = cputime;
			[C,P,Q,score] = self.greedy(G);
			Qs = Q;
			dt = cputime-dt;
		end
		function [q,qs,score] = eval(self,A,x,varargin)
			q = 1-(sum(x)-max(x))*2/max(1,length(x)-2);
			qs = q;score =x;
		end
		function param = initParam(self,param)
			if ~isfield(param,'disp');param.disp = false;end
		end
	end
	methods (Access = private)
		function [C,P,Q,score] = greedy(self,G,param)
			n = G.numNode();
			d = G.degree();
			m = G.numEdge();
			A = G.adjacency_matrix('binary');
			
			
			x = zeros(n,1);
			[idx,val] = min(d);
			r = find( d==idx);
			idx = randsample(1:length(r),1);
			idx = r(idx);x(idx) = 1;
			ak = d(idx);bk = 0;
			alpha = zeros(n,1);
			alpha(1) = 0;
			
			for k = 2:n
				score = (2*ak * (x'*A)'  - bk * d)./max(1,ak*(ak+d));	
				score(x>0) = Inf;	
				idx = min(score);	
				r = find( score ==idx );
				idx = randsample(1:length(r),1);
				idx = r(idx);
				x(idx) = 1;
				ak = ak + d(idx);bk = x'*A*x;
				alpha(idx) = bk/max(1,ak);
			end
			C = alpha;
			P = 1-alpha;
			Q = eval(self,A,alpha);	
			score = alpha;
		end
	end
end
