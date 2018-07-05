classdef eigen_cp 
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
			g = graph(A);
			C = abs(g.max_eigen());
			P = 1-C;
			score = C;Q = C'*A*C;Qs = Q;
		end
		function q = eval(self,A,x)
			q = x'*A*x;
		end
	end
	methods (Access = private)
		function model = init_model(self,model)
			if ~isfield(model,'solver');model.solver = 'eigen';end
		end
		
	end
end
