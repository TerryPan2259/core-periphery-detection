classdef fabio_cp 
	methods (Access = public)
		function [C,P,Q,score] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
		
			switch model.solver
				case 'greedy'
					[C,P,Q,score] = self.greedy(graph(A),model);
				otherwise	
					disp('unknown solver for borgatt & Everett model')
			end
		end
	end
	methods (Access = private)
	
		function [C,P,Q,score] = greedy(self,G,model)
			cmp = G.connected_components();
			N = G.N;
			% compute minimum eigen vector of modularity matrix for each components 
			Q = [];
			C = [];	
			P = [];	
			for c = cmp 
				At= G.A(c>0,c>0);
				sa = G.deg(c>0); 
				% ----------------------------	
				% optimise	
				% ----------------------------
				Cs = [];
				for nid = 1:size(At,1) 
					Cn = sparse(nid,1,1,size(At,1),1);
					while(true)
						score = (At *Cn);
						score(Cn>0)  = intmax;
						[val,idx] = min(score);
						Cn(idx) = 1;
						Q = Cn'*At*Cn / sum(sa(Cn>0));
						if Q > model.alpha
							Cn(idx) = 0;
							break
						end	
					end
				
					if sum(Cn(:)) >= sum(Cs(:))
						Cs = Cn;
					end
				end
				% ----------------------------	
				tmp = sparse(G.N,1,0);
				tmp(c>0) = Cs;	
				Qt = Cs'*At*Cs / sum(sa(Cs>0));
				C = [C,tmp];
				
				tmp = sparse(N,1,0);
				tmp(c>0) = 1-Cs;
				P = [P,tmp];
				
				Q =[Q , Qt];	
			end
			C = sum(C,2);
			P = 1-C;
			if size(C,2)==0
				Q = 0;
				C = sparse(G.sz(),1,0);
				P = 1-C;score = C;
				return;
			end
			
			score = ones(G.sz(),1);
		end
		
		function model = init_model(self,model)
			if ~isfield(model,'alpha');model.alpha = 0.1;end
			if ~isfield(model,'solver');model.solver = 'greedy';end
		end
	end
end
