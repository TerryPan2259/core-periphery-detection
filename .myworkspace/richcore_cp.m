classdef richcore_cp < cpabst
	methods (Access = public)
		function [C,P,Q,Qs,score,param,dt] = detect(self, G, param)
			A = G.adjacency_matrix();
			C= [];  P = [];Q = [];score = [];
			dt = cputime;	
			deg = sum(A,2);
			n = size(A,1);
			[~,ord] = sort(deg,'descend');
			x = sparse(n,1,0);
			kpmax = -Inf; 
			for i = 1:n
				nid = ord(i);
				kp = A(nid,:)*x;	
				x(nid) = 1;	
				if kpmax <kp
					C = x;
					kpmax = kp;
				end
			end
			P = 1-C;
			Q = self.eval(G,C);Qs = Q;
			dt = cputime-dt;	
		end
				
		function [q,qs,score] = eval(self,G,x,varargin)
			A = G.adjacency_matrix();
			deg = G.degree();
			q = x'*A*x/sum(deg.*x);
			qs = q;score =x;
		end
		
		function param = initParam(self,param)
			if isempty( param )	
				param.name = 'richcore';
			end	
			if ~isfield( param, 'disp' ); param.disp = false; end
		end
	end
end	
