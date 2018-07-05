classdef richcore
	methods (Access = public)
		function [C,P,Q,Qs,score,model] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			C= [];  P = [];Q = [];score = [];
			
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
			Q = self.eval(A,C);Qs = Q;
		end
				
		function q = eval(self,A,x)
			deg = sum(A,2);
			q = x'*A*x/sum(deg.*x);
		end
	end
	methods (Access = private)
		
		
		function model = init_model(self,model)
		end
		
	end
end	
