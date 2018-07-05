% this depends on community detection methods
classdef modular2_cp 
	methods (Access = public)
		function [C,P,Q,score] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
			
			com = community();
			
			% community detection step
			mymodel.name = 'modularity';
			U = com.detect(A,mymodel);
			
			% periphery detection step
			com = community();
			mymodel.name = 'modularity_conductance';
			C = [];P = [];
			n = size(A,1);	
			for i = 1:size(U,2)
				u = U(:,i);
				sz = sum(u);
				As = A(u>0,u>0);
				R = sparse(find(u),1:sz,1,n,sz);
			
				D = com.detect(As,mymodel);
				C = [C,R*D(:,1)]; P = [P,R*D(:,2)]; 
			end
			Q = 0;score = zeros(n,1);Qs = zeros(1,size(C,2));
		end
	end
	methods (Access = private)
		function model = init_model(self,model)
			if ~isfield(model,'itnum');model.itnum = 10;end
		end
	end
end
