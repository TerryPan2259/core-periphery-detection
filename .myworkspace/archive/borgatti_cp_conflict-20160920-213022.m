classdef borgatti_cp 
	methods (Access = public)
		function [C,P,Q,score] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
			C= [];  P = [];Q = [];score = [];
			switch model.solver
				case 'eigen'
					[C,P,Q,score] = self.eigen(graph(A));
				otherwise	
					disp('unknown solver for borgatt & Everett model')
			end
		end
	end
	methods (Access = private)
	
		function [C,P,Q,score] = eigen(self,G)
			cmp = G.connected_components();
			N = G.N;
			% compute minimum eigen vector of modularity matrix for each components 
			Q = [];
			C = [];P=[];	
			for c = cmp 
				At= G.A(c>0,c>0);
				v = G.maximum_eigen_vector(At);
				v = abs(v);
				[vals,idx] = sort(v,'descend');
				
				xbest = sparse(length(v),1,0);
				bestscore = 0;	
				for i = 1:length(v)
					xcand = sparse(length(v),1,0);xcand(v>=vals(i))=1;
					template = xcand*xcand'+ xcand*(1-xcand)' + (1-xcand)*xcand';
					template = template - diag(diag(template));
					score = corr(At(:),template(:));
					if score > bestscore
						bestscore = score;
						xbest = xcand;
					end	
				end	
				tmp = sparse(N,1,0);
				tmp(c>0) = xbest;
					
				C = [C,tmp];
				
				tmp = sparse(N,1,0);
				tmp(c>0) = 1-xbest;
				P = [P,tmp];
				Q = [Q, bestscore];	
			end
			C = sum(C,2);
			P = 1-C;
			
			if size(C,2)==0
				Q = 0;
				C = sparse(G.sz(),1,0);
				P = 1-C;score = C;
				return;
			end
			score = sum(G.A,2);
		end
		
		function model = init_model(self,model)
			if ~isfield(model,'solver');model.solver = 'eigen';end
		end
		
	end
end
