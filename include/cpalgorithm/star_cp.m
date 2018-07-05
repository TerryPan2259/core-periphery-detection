classdef star_cp < cpabst
	
	methods ( Access = public )
		
		function [C, P, Q, Qs, score, param, cpu_time] = detect( self, G, varargin )
			% --------------------------------
			% Initialise
			% --------------------------------
			param = [];
			if nargin == 3
				param = varargin{1};
			end
			param = self.initParam( param );
			C = []; P = []; Qs = []; score = []; Q = - Inf;
			cpu_time = cputime;
	
			% --------------------------------
			% Optimise node labels using a label switching algorithm 
			% --------------------------------
			for it = 1:param.numRun 
				if param.overlap 	
					if strcmp(param.optalgorithm,'greedy') 
						[Ct, Pt, Qt, Qst] = self.greedy_overlap( G );
					elseif strcmp(param.optalgorithm,'ip')
						[Ct, Pt, Qt, Qst] = self.ip_overlap( G );
					end
				else
					if strcmp(param.optalgorithm,'greedy') 
						[Ct, Pt, Qt, Qst] = self.greedy_exclusive( G );
					elseif strcmp(param.optalgorithm,'ip')
						[Ct, Pt, Qt, Qst] = self.ip_exclusive( G );
					end
				end
				
				if Qt > Q
					C = Ct; P = Pt; Q = Qt; Qs = Qst; 
				end
			end
			cpu_time = cputime - cpu_time;
		end
		
		function [Q, Qs, score] = eval( self, G, C, P, varargin )
			Q = size(C,2);
			Qs = ones(1,size(C,2));
			score =sum(C,2);
			score = Q*score/sum(score); 
		end
	
		function param = initParam( self, param )
			if isempty( param )	
				param.name = 'star';
			end	
			if ~isfield( param, 'numRun' ); param.numRun = 20; end
			if ~isfield( param, 'overlap' ); param.overlap = false; end
			if ~isfield( param, 'optalgorithm' ); param.optalgorithm = 'greedy'; end
			if ~isfield( param, 'disp' ); param.disp = false; end
		end
	end
	
	methods ( Access = private )
		
		% Aslam, J., Pelekhov, E., Rus, D.  J. Graph. Algo. Appl. 2004 
		function [C,P,Q,Qs,score] = greedy_overlap(self,G)
			N = G.numNode();
			deg = G.degree();
			remain = ones(N,1);
			A = G.adjacency_matrix('binary');
		
			C = [];P = [];
			while( nnz(remain) ~= 0 )
				[maxdeg,idx] = max(remain.*deg);
			
				c = sparse(idx,1,1,N,1);
				p = sparse(N,1,0);p( remain==1 & A(:,idx)==1 ) = 1;
				remain(idx) = 0;
				remain( p==1 ) = 0;
				
				C = [C,c];P = [P,p];
			end
			
			Q = size(C,2);
			Qs = ones(1,size(C,2));
			score =sum(C,2);
			score = Q*score/sum(score); 
		end
		function [C,P,Q,Qs,score] = ip_overlap(self,G)
			N = G.numNode();
			deg = G.degree();
			remain = ones(N,1);
			A = G.adjacency_matrix('binary');
			
			[r,c] = find(triu(A));
			x = binvar(2*length(r),1);
			X = sparse([r;c],[c;r],x);
			y = binvar(N,1);
			
			Const = [sum(X,2) + y >= 1, x<=[y(c);y(r)]];
			Obj = sum(y);
		
			options = sdpsettings('verbose',0);
				
			sol = optimize(Const, Obj, options);
			x = value(x);
			y = value(y);
				
			X = sparse([r;c],[c;r],x);
			centre = find(any(X,1));
			P = X(:,centre);
			C = sparse(centre,1:length(centre),1,N,length(centre));
			
			
			Q = size(C,2);
			Qs = ones(1,size(C,2));
			score =sum(C,2);
			score = Q*score/sum(score);
		end
		
		function [C,P,Q,Qs,score] = greedy_exclusive(self,G)
			N = G.numNode();
			deg = G.degree();
			remain = ones(N,1);
			A = G.adjacency_matrix('binary');
		
			C = [];P = [];
			while( nnz(remain) ~= 0 )
				[maxdeg,idx] = max(deg);
				
				if maxdeg ==0
					break;
				end
			
				nei = find(A(:,idx)>0);

			
				c = sparse(idx,1,1,N,1);
				p = sparse(N,1,0);p( nei ) = 1;
				C = [C,c];P = [P,p];
				
				
				remain( idx ) = 0;
				remain( nei ) = 0;
				
				A([idx;nei],:) = 0;	
				A(:,[idx;nei]) = 0;	
				
				deg = sum(A,2);
			end
			
			remain = find(remain);
			for i = remain'
				C(i,size(C,2)+1) = 1;
				P(N,size(P,2)+1) = 0;
			end
			
			Q = size(C,2);
			Qs = ones(1,size(C,2));
			score =sum(C,2);
			score = Q*score/sum(score); 
		end
		
		function [C,P,Q,Qs,score] = ip_exclusive(self,G)
			N = G.numNode();
			deg = G.degree();
			remain = ones(N,1);
			A = G.adjacency_matrix('binary');
			
			[r,c] = find(triu(A));
			x = binvar(2*length(r),1);
			X = sparse([r;c],[c;r],x);
			y = binvar(N,1);
			
			Const = [sum(X,2) + y == 1, x<=[y(c);y(r)]];
			Obj = sum(y);
		
			options = sdpsettings('verbose',0);
				
			sol = optimize(Const, Obj, options)
			x = value(x);
			y = value(y);
				
			X = sparse([r;c],[c;r],x);
			centre = find(any(X,1));
			P = X(:,centre);
			C = sparse(centre,1:length(centre),1,N,length(centre));
			
			
			Q = size(C,2);
			Qs = ones(1,size(C,2));
			score =sum(C,2);
			score = Q*score/sum(score); 
		end
	end
end
