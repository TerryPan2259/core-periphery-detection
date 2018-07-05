% Two-step algorithm
%
% @version 0.1.0
% @date 2017/02/19
%
% @author Sadamori Kojaku 
classdef two_step_cp 

	methods ( Access = public )

		function [C, P, Q, Qs, score, param, cpu_time] = detect( self, G, varargin )
			param = [];
			if nargin == 3
				param = varargin{1};
			end
			if strcmp(param.name,'two-step')
				[C, P, Q, Qs, score, param, cpu_time] = self.detect_individually(G,param);	
			elseif strcmp(param.name,'divisive')
				[C, P, Q, Qs, score, param, cpu_time] = self.detect_recursively(G,param); % bad name really...	
			end
		end
		
		function [q, qs, score] = eval( self, G, C, P, varargin )
			deg = G.degree();
			M = G.numEdge();
			A = G.adjacency_matrix();
			U = C+P;
			degU = sum( diag(deg) * U , 1 );
			
			qs = ( sum( (A*U).*U, 1 ) - degU.*degU/(2*M) )/(2*M);	
			q  = sum(qs);
			score = (q/size(A,1)) * ones(size(A,1),1);	
		end
		
		function param = initParam( self, param )
			if ~isfield( param, 'name' ); param.name='two-step'; end
			if ~isfield( param, 'disp' ); param.disp = false; end
			if ~isfield( param, 'numRun' ); param.numRun = 20; end
		end
	end
	
	methods ( Access = private )
		% detect communities and core-periphery individually, then merge the results 
		function [C, P, Q, Qs, score, param, cpu_time] = detect_individually( self, G, param )
			% --------------------------------
			% Initialise
			% --------------------------------
			param = self.initParam( param );
			C= []; P = []; Q = -Inf; score = [];
			cpu_time = cputime;	
			
			% --------------------------------
			% Detect a single core-periphery pair and communities  
			% --------------------------------
			cpd = cpalgorithm();
			beparam = cpd.init('be');
			beparam.numRun = param.numRun;
			[C, P] = cpd.detect( G.adjacency_matrix(), beparam );
			
			% detect communities by modularity	
			cm = comalgorithm();
			cmparam = cm.init('modularity');
			cmparam.numRun = param.numRun;
			[U,Q,Qs,score] = cm.detect( G.adjacency_matrix(), cmparam );
			N =G.numNode();
			if isempty(U)
				P = sparse(N,1,0);
				C = sparse(N,1,0);
				Q = 0;Qs = 0;score = sparse(N,1,0);
				cpu_time = cputime - cpu_time;
				
				return
			end	
			C = sparse( 1:size( C, 1 ), 1:size( C, 1 ), C ) * U;
			P = U - C;
			%[Q, Qs, score] = self.eval( G, C, P );
			cpu_time = cputime - cpu_time;
				
		end
		
		% detect communities first, then divide each of them into core and periphery
		function [C, P, Q, Qs, score, param, cpu_time] = detect_recursively( self, G, param )
			% --------------------------------
			% Initialise
			% --------------------------------
			param = self.initParam( param );
			C= []; P = []; Q = -Inf; score = [];
			cpu_time = cputime;	
			A = G.adjacency_matrix();
			
			% --------------------------------
			% detect communities by modularity	
			% --------------------------------
			cm = comalgorithm();
			cmparam = cm.init('modularity');
			cmparam.numRun = param.numRun;
			[U,Q,Qs,score] = cm.detect( A, cmparam );
			N =G.numNode();
			if isempty(U)
				P = sparse(N,1,0);
				C = sparse(N,1,0);
				Q = 0;Qs = 0;score = sparse(N,1,0);
				cpu_time = cputime - cpu_time;
				
				return
			end	
			
			% --------------------------------
			% detect cp in each community 
			% --------------------------------
			cpd = cpalgorithm();
			beparam = cpd.init('be');
			beparam.numRun = param.numRun;
			N = size(A,1);
			C = [];P = [];
			for i = 1:size(U,2)
				slice = find(U(:,i));
				R = sparse(slice,1:length(slice),1,N,length(slice));
				B = A(slice,slice);
				[Ci, Pi] = cpd.detect( B, beparam );
				if isempty(Ci) | length(Ci)==0
					
					C = [C,U(:,i)];
					P(:,size(P,2)+1) = 0;
				else
					C = [C,R*Ci];
					P = [P,R*Pi];
				end
			end
			cpu_time = cputime - cpu_time;
				
		end
	end

end
