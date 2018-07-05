classdef borgatti_multi
	methods (Access = public)
		function [C,P,Q,Qs,score,model] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			
			model = self.init_model(model);
			cpd = borgatti_cp();	
			[C,P,Q,Qs,score] = cpd.detect(A,model);	
				
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
		end
	end
	
end
