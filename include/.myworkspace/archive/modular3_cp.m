% this depends on community detection methods
classdef modular3_cp 
	methods (Access = public)
		function [C,P,Q,score] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
	
			% degree cut -------------------
			G = graph(A);
			d = G.degree();;
			m = sum(d)/2;
			dave = 2*m/G.sz();
			below = d < dave;
			Ag = A;
			Ag(below,below) = 0;
			M = graph(Ag).modularity_matrix(sum(Ag(:))/sum(A(:)));
			%M(below,below) = -d(below)*d(below)'/(4*m*m);
			
			% community decomposition ----- 
			com = community();
			mymodel.name = 'gpotts';
			mymodel.solver = 'louvain';
			mymodel.M = M;
			U = com.detect(Ag,mymodel);
			
			isolated = sum(A,1)==0;
			 
			%U(below,:) = 0;U=U(:,any(U,1));
			
			% core/periphrey decomposition -
			com = cpdetection();
			mymodel.name = 'borgatti';
			mymodel.solver = 'eigen';
			C = [];P = [];
			n = size(A,1);	
			for i = 1:size(U,2)
				u = U(:,i);
				sz = sum(u);
				As = Ag(u>0,u>0);
				R = sparse(find(u),1:sz,1,n,sz);
				
				[D,dQ]= com.detect(As,mymodel);
				if length(D)==0
					C = [C,R*ones(size(As,1),1)]; P = [P,zeros(size(A,1),1)]; 
				else
					C = [C,R*D]; P = [P,R*(1-D)]; 
				end
			end
	
			% periphery merge ------------
			rest = find(isolated);
			Mc = M*C;
			for nid = rest
				[dq,cid] = max(Mc(nid,:));
				if dq <=eps*1000;continue;end
				P(nid,cid) = 1;
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
