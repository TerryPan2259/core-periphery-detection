classdef xian_cp < cpabst

	methods ( Access = public )

		function [C, P, Q, Qs, score, param, cpu_time] = detect( self, G, param )
			% --------------------------------
			% Initialise
			% --------------------------------
			C= []; P = []; Q = -Inf; score = [];Qs = [];
			ts = cputime;
			
			% --------------------------------
			% Optimise node labels using the Kernighan-Lin algorithm 
			% --------------------------------
			[C,P,Q,Qs,score, param] = self.detect_multiple_cp(G,param);
		
			if param.exclusive
				iscore = any(C,2);
				overlap = find( sum(C+P,2)> 1 )';
				for nid = overlap 
					if iscore(nid)
						cid = randsample(size(C,2),1,true,full(C(nid,:)))
						C(nid,:) = 0;C(nid,cid) = 1;	
					else
						cid = randsample(size(P,2),1,true,full(P(nid,:)))
						P(nid,:) = 0;P(nid,cid) = 1;	
					end	
				end
			end
			
			cpu_time = cputime-ts;
		end
		
		function [q,qs,score] = eval(self,G,x,varargin)
			q = NaN;qs = NaN * ones(1,size(x,2));score = NaN*ones(size(x,1),1);
		end
			
		function param = initParam(self,param)
			if ~isfield(param,'numRun');param.numRun = 20;end
			if ~isfield(param,'name') param.name = 'be';end
			if ~isfield(param,'beta') param.beta = 1-eps*1000;end
			if ~isfield(param,'disp') param.disp = false;end
			if ~isfield(param,'exclusive') param.exclusive = true;end
		end
	end
	
	methods ( Access = private )

		function [C, P, Q, Qs, score, param] = detect_multiple_cp( self, G, param )
			
			A = G.adjacency_matrix();
			N = size(A,1);
			deg = G.degree();
			
			% fint an initial node
			
			%clos=zeros(N,1);  % initialize closeness vector
			%[r,c] = find(triu(A,1));
			%D = dijkstra([r,c],size(A,1),length(r));
			%for i=1:N; clos(i)=1/sum( D(i,:) );end
			%[~,nid] = max(clos);
			[maxdeg,nid] = max(deg);	
			U = sparse(N,1,0);U(nid) = 1;
			nodeorder = zeros(N,1);nodeorder(1) = nid;
			
			ordidx = 2;	
			while nnz(U) < N
					
				score = A*U + deg/maxdeg;
				score(U>0) = 0;
				
				[val,idx] = max(score);
				nid = randsample( N, 1, true, full(double(val == score)) );
				
				U(nid) = 1;	
				nodeorder(ordidx) = nid;
				ordidx = ordidx + 1;	
			end
			C = [];P =[]; Q = [];Qs = [];
			
			alpha = floor(mean(deg));RD = zeros(N,1);
			for i = 1:N
				
				ns = 1;
				if alpha < i
					ns = i-alpha+1;
				end
				 
				As = A( nodeorder(ns:i), nodeorder(ns:i) );
				
				RD(i) = nnz(As)/(size(As,1)*(size(As,1)-1));	
			end
			
			RD(isnan(RD)) = 0;
		
			% find core set
			C = [];c = sparse(N,1);
			for i = alpha:N
				if RD(i) >= param.beta 
					c( nodeorder( (i-alpha + 1):i) ) = 1;	
				else
					if (RD(i-1) >= param.beta - 1e-5) & i > alpha
						C = [C,c];c = sparse(N,1);
					end	
				end	
			end
			if RD(N) >= param.beta
				C = [C,c];
			end
		
			param.RD = RD;
			if isempty(C) % if no core is found
				C = sparse(N,1,0);P = sparse(N,1,0);
				P(1:N) = 1;
				score = NaN;
				Q = NaN;Qs = NaN*ones(1,size(C,2));
				return;
			end
			
			% find periphery
			P = sparse(N,size(C,2),0);	
			while true
				V = find(~any(C+P,2))';
				if isempty(V);break;end
				for nid = V
					score = A(nid,:) * (C+P);	
					[val,cid] = max(score);
					P(nid,score==val) = 1;	
				end
			end
		
			% find core-periphery set	
			iscore = any(C,2);
			overlap = sum(C,2)>1 | sum(P,2) > 1;
			C(overlap,:) = 0; P(overlap,:) = 0;
			V = find(overlap)';
			for nid = V
				score = A(nid,:) * (C+P);	
				[val,cid] = max(score);
				if iscore(nid)
					C(nid,val==score) = 1;	
				else
					P(nid,val==score) = 1;	
				end
			end
			score = NaN;
			Q = NaN;Qs = NaN*ones(1,size(C,2));
			param.RD = RD;
		end
	end
end
