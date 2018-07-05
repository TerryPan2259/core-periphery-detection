classdef rombach_multi
	methods (Access = public)
		function [C,P,Q,Qs,score,model] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			
			model = self.init_model(model);
			cpd = rombach_cp();	
			[C,P,Q,score] = cpd.detect(A,model);	
			
			cmd = community();	
		
			clear cmodel;
			cmodel.name = 'modularity';	
			[U,~,Qs] = cmd.detect(A,cmodel);
			C = sparse(1:size(C,1),1:size(C,1),C)*U;
			P = U-C;
		end
	end
	methods (Access = private)
					
		function model = init_model(self,model)
			if ~isfield(model,'disp');model.disp = false;end
			if ~isfield(model,'alpha');model.alpha = 1;end
			if ~isfield(model,'beta');model.beta = 0.5;end
			%if ~isfield(model,'solver');model.solver = 'label_switching';end
			if ~isfield(model,'solver');model.solver = 'annealing';end
		end
	end
	
end
