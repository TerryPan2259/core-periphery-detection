classdef rombach_cp 
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
				case 'eigen'
				case 'label_switching'
					[C,P,Q,score] = self.label_switching(A,model.alpha,model.beta);
				case 'annealing'
					[C,P,Q,score] = self.simulated_annealing(A,model.alpha,model.beta);
				otherwise	
					disp('unknown solver for Rombach method')
			end
			Qs = Q;
		end
	end
	methods (Access = private)
		
		function [C,P,Q,score] = simulated_annealing(self,A,alpha,beta)
			n = size(A,1);
			x = randsample(n,n);
			c = self.corevector(x,alpha,beta);	
				
			function x = gen(x)
				nids = randsample(n,2);
				x(nids) = x(nids(2:-1:1));	
			end	
			clear aopt
			aopt.Generator = @gen;
			aopt.Verbosity = 0;
			[C,fval] = anneal(@(x)(-self.eval(A,x)),c,aopt);
			P = 1-C;
			score = C;
			Q = -fval;
		
		end
		
		function q = eval(self,A,c)
			q = c'*A*c; 
		end
		
		function c = corevector(self,x,alpha,beta)
			n = length(x);
			bn = floor(beta*n);
		
			cx = x<=bn;
			px = ~cx;
			
			c = (1-alpha)/(2*bn)*x.*cx + ((x.*px-bn)*(1-alpha)/(2*(n-bn))  + (1+alpha)/2).*px;
		end
					
		function model = init_model(self,model)
			if ~isfield(model,'disp');model.disp = false;end
			if ~isfield(model,'alpha');model.alpha = 0.8;end
			if ~isfield(model,'beta');model.beta = 0.5;end
			%if ~isfield(model,'solver');model.solver = 'label_switching';end
			if ~isfield(model,'solver');model.solver = 'annealing';end
		end
	end
end
