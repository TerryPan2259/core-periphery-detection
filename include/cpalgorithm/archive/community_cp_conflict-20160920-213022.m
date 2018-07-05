% this class depends on /program/lib/community
classdef community_cp 
	methods (Access = public)
		 function [C,P,Q,score] = detect(varargin)
			% -----------------------------------
			% Initilise 
			% -----------------------------------
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
			
			switch model.cpmodel.name
				case 'borgatti'
					cpd = borgatti_cp();
				case 'fabio'
					cpd = fabio_cp();	
				otherwise
					display ('unknown model name in recursive_cp')
			end
			
			% -----------------------------------
			% community division
			% -----------------------------------
			com = community();
			[Com,Q] = com.detect(A,model.communitymodel);
			
			if size(Com,1)==0
				score  = sparse(1:size(A,1),1,0);
				C = sparse(1:size(A,1),1,1);P=1-C;
				return
			end	
			% -----------------------------------
			% core/periphery division
			% -----------------------------------
			C = [];P=[];
			for k = 1:size(Com,2)
				At = A;
				At(Com(:,k)==0,:) = 0; 
				At(:,Com(:,k)==0) = 0; 
				
				[Ct,~,Qt] = cpd.detect(At,model.cpmodel);
				
				
				if size(Ct,1)==0;C = [C,Com(:,k)];P = [P,sparse(size(C,1),1,0)];break;end
				C = [C,Ct];
				Pt = Com(:,k).*(1-Ct);
				P = [P,Pt];
				Q = [Q,Qt];
			end
			
			score = ones(size(C,1),1);
		end
	end
	
	methods (Access = private)
		function model = init_model(self,model)
			if ~isfield(model,'cpmodel');model.cpmodel.name = 'borgatti';end
			if ~isfield(model,'communitymodel');model.communitymodel.name = 'modularity';end
		end
	end
end	
