classdef recursive_cp 
	methods (Access = public)
		 function [C,P,Q,score] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
			
			switch model.recursivemodel.name
				case 'borgatti'
					cpd = borgatti_cp();
				case 'fabio'
					cpd = fabio_cp();	
				otherwise
					display ('unknown model name in recursive_cp')
			end
	
			remain = true(size(A,1),1);
			At = A;C = [];	Q = [];
			for k = 1:model.K	
				[Ct,~,Qt] = cpd.detect(At,model.recursivemodel);
				if size(Ct,1)==0;break;end
				remain = remain & Ct(:,1)==0;
				C = [C,Ct(:,1)];
				At(~remain,:) = 0;
				At(:,~remain) = 0;
				
				if any(remain)==false
					break;
				end
				Q = [Q,Qt];
			end
			
			score = ones(size(C,1),1);
			
			% find periphery	
			P = A*C;% - d*sum(C,1)/(2*m);
			P(P<0)=0;	
			if size(C,2)==1
				P =sign(P); 
				P(any(C,2)) =0;
				return
			end
				
			[val,cid] = max(P');
			
			P = sparse(1:size(A,1),cid,1,size(A,1),size(C,2));
			P(any(C,2) | val'==0,:) =0;
			
		end
	end
	
	methods (Access = private)
		function model = init_model(self,model)
			if ~isfield(model,'recursivemodel');model.recursivemodel.name = 'borgatti';end
			if ~isfield(model,'K');model.K = 2;end
		end
	end
end	
