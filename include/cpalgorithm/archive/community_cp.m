% this class depends on /program/lib/community
classdef community_cp 
	methods (Access = public)
		 function [C,P,Q,Qs,score,model] = detect(varargin)
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
						
			% -----------------------------------
			% community division
			% -----------------------------------
			com = community();
			model.communitymodel
			[Com,Q,Qs,~] = com.detect(A,model.communitymodel);
			if model.disp	
				disp(sprintf('%d communities detected',size(Com,2)));
			end	
			if size(Com,1)==0
				C = sparse(1:n,1,1);P=1-C;
				return
			end	
			% -----------------------------------
			% core/periphery division
			% -----------------------------------
			isolated = sum(Com,1)<=1;
			Com(:,isolated) = [];
			Qs(:,isolated) = [];
			cpd = cpdetection();
			C = [];P=[];
			for k = 1:size(Com,2)
				slice = Com(:,k)>0;
				At = A(slice,slice);
				
				[Ct,~,Qt] = cpd.detect(At,model.cpmodel);
				
				c = zeros(size(Com,1),1);
				p = zeros(size(Com,1),1);
				if size(Ct,1)~=1
					c(slice) = Ct;
					p(slice) = 1-Ct;
				else
					p(slice)=1;
				end
					
				C = [C,c];
				P = [P,p];
				Q = [Q,Qt];
				if model.disp	
					disp(sprintf('%d/%d core/periphery detected',k,size(Com,2)));
				end	
			end
			isolated = find(~any(Com,2))';
			if length(isolated)>0
				P(:,size(C,2)+1) = 0;
				P(isolated,size(C,2)+1) = 1;
				Qs(size(C,2)+1) = 0;	
				C(:,size(C,2)+1) = 0;	
			end
			if size(C,1)~=size(A,1);
				C = zeros(size(A,1),1);
				P = zeros(size(A,1),1);
				Qs = zeros(1,size(C,2));	
				score = ones(size(C,1),1);
				return
			end
			score = ones(size(C,1),1);
		end
	end
	
	methods (Access = private)
		function model = init_model(self,model)
			if ~isfield(model,'cpmodel');model.cpmodel.name = 'borgatti';end
			if ~isfield(model,'disp');model.disp = false;end
			if ~isfield(model,'communitymodel');model.communitymodel.name = 'modularity';end
		end
	end
end	
