classdef normalised_cp 
	methods (Access = public)
		function [C,P,Q,score] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
			
			d = sum(A,2);
			m = sum(d)/2;
			rho = 1;
			
			while(true)
				M = (A - rho * d*d'/(2*m))/(2*m);
				C = self.eigen(M);
				C
				Q = C'*M*(1-C);
				if Q ==0
					break;
				end
				rho = 2*m * C'*A*(1-C)/( C'*d*((1-C)'*d) );
			end
			P = 1-C;
			Q= 2*m * C'*A*(1-C)/( C'*d*((1-C)'*d) );
			score = 0;
		end
	end
	methods (Access = private)
	
		function C = eigen(self,M)
			N = size(M,1);
			% calc stationary distribution and connected components
			C = sign(graph(M).minimum_eigen_vector(M))>0;
		end
		function C = opt(self,M)
			N = size(M,1);
			% calc stationary distribution and connected components
			X = sdpvar(N,1);
			th = sdpvar(1);
			F = [th >= trace(X'*M*X)];
			F = [F,binary(X)];
			optimize(F,th);
			C = value(X);
		end
		
		function model = init_model(self,model)
			if ~isfield(model,'alpha');model.alpha = 0.1;end
			if ~isfield(model,'solver');model.solver = 'greedy';end
		end
	end
end
