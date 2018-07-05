classdef yan_cp < cpabst

	methods ( Access = public )

		function [C, P, Q, Qs, score, param, cpu_time] = detect( self, G, param )
			% Initialise ----------
			C= []; P = []; Q = -Inf; score = [];Qs = [];
			A = G.adjacency_matrix('binary');
			N = size(A,1);
			ts = cputime;
			
			% compute distance between nodes
			D = dijkstra(A,(1:N)');	
			% find hierarchical clusters using average linkage algorithm
			H = linkage(D(triu(true(N),1))','average');
			U = sparse(1:N,1:N,1); % store all clusters
			cids = sparse(2*N,1,0);
			cids(1:N) = 1;
			s = {};  % store indexes of clusters at each aggromerative level
			for i = 1:size(H,1)
				sid = H(i,1);
				did = H(i,2);
				newcid = size(U,2) + 1;
				U(:,newcid) = U(:,sid) + U(:,did);
				cids(sid) = 0;	
				cids(did) = 0;
				cids(newcid) = 1;
				s{i} = find(cids)';	
			end	
			
			% compute r-score (core-periphery-ness)
			r = [];
			z = [];
			deg = sum(A,2);
			M = sum(deg)/2;
			for i = 1:size(H,1)
				Unow = U(:,s{i});
				if(size(Unow,2)==1); continue;end;
					
				% divide them into cores and a periphery
				iscore = sum(Unow,1)>1;
				C = Unow(:,iscore);
				P = double(any(Unow(:,~iscore),2));	
					
				% density of cores
				r(i,1) = self.eval( G, C, P, param);
				
				if r(i,1)<eps
					z(i,1) = 0; 
				else
					Wcc = trace(C'*A*C);
					EWcc = sum((deg'*C).^2)/(2*M);
					sigWcc = sqrt(EWcc);	
				
					z(i,1) = (Wcc-EWcc) / sigWcc; 
				end
				
			end
			[r,z]
			[val,idx] = max(z);
				
			Unow = U(:,s{idx});
			iscore = sum(Unow,1)>1;
			C = Unow(:,iscore);
			P = double(any(Unow(:,~iscore),2));
			tmp = sparse(N,size(C,2),0);
			tmp(:,1) = P;
			P = tmp;
				
			Q = r(idx);
			Qs = Q * ones(1,size(C,2));
			score = Q  * ones(size(C,1),1);
			
			cpu_time = cputime-ts;
		end
		
		function [q,qs,score] = eval(self,G,C,P,param)
			A = G.adjacency_matrix('binary');
			N = size(A,1);
			M = sum(A(:))/2;
				
			n = sum(C,1);
			Wcc = trace(C'*A*C);
			barWcc = 2*sum(sum(A*P,2)) - trace(P'*A*P);
			dcore =  Wcc / sum(n.*(n-1));
			dperi = barWcc /(N*(N-1) - sum(n.*(n-1))); 
			
			q = dcore / dperi;
			if isnan(q) | isinf(q)
				q = 0;
			end
				
			qs = q * ones(1,size(C,2))/size(C,2);score = q*ones(size(N,1),1)/N;
		end
			
		function param = initParam(self,param)
			addpath /panfs/panasas01/emat/sk16722/program/lib/dijkstra
			addpath /panfs/panasas01/emat/sk16722/program/lib/randgraph
			
			if ~isfield(param,'numRun');param.numRun = 1;end
			if ~isfield(param,'name') param.name = 'yan';end
			if ~isfield(param,'disp') param.disp = false;end
		end
	end
	
	methods ( Access = private )

	end
end
