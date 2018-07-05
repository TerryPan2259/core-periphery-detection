% Discrete variant of the Borgatti-Everett algorithm 
%
% @version 0.1.0
% @date 2017/02/19
%
% @author Sadamori Kojaku 
classdef be_cp < cpabst

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
			if param.dcor  % degree correlation
				for it = 1:param.numRun
					if strcmp(param.name,'be_cont_dcor') | strcmp(param.name,'be_continuous')
						[Ct,Pt,Qt,Qst,scoret] = self.becont_dcor2(G,param);
					%	[Ct,Pt,Qt,Qst,scoret] = self.becont_dcor(G,param);
					end
							
					if Qt > Q | it == 1
						C = Ct; P = Pt; Q = Qt; score = scoret; Qs = Qst;
					end
				end
			else
				for it = 1:param.numRun
					if strcmp(param.name,'be') 
						[Ct,Pt,Qt,Qst,scoret] = self.Kernighan_Lin_algorithm_mex(G,param);
						%[Ct,Pt,Qt,Qst,scoret] = self.Kernighan_Lin_algorithm_mex(G,param);
					elseif strcmp(param.name,'be_without_cp_edges')
						[Ct,Pt,Qt,Qst,scoret] = self.Kernighan_Lin_algorithm_without_cp_edges_mex(G,param);
					elseif strcmp(param.name,'be_continuous')
						[Ct,Pt,Qt,Qst,scoret] = self.becontinuous(G,param);
					end
							
					if Qt > Q | it == 1
						C = Ct; P = Pt; Q = Qt; score = scoret; Qs = Qst;
					end
				end
			end
			cpu_time = cputime-ts;
		end
		
		function [q,qs,score] = eval(self,G,x,varargin)
			name = 'be';
			dcor = false;	
			if nargin ==5
				param = varargin{2};
				if isfield(param,'name') name = param.name;end
				if isfield(param,'dcor') dcor = param.dcor;end
			end
			
			if dcor 
				if strcmp(name,'be_cont_dcor')
					[q,qs,score] = self.eval_cont_dcor(G,x);	
				end
			else
				if strcmp(name,'be')
					[q,qs,score] = self.eval_with_cp(G,x);	
				elseif strcmp(name,'be_without_cp_edges')
					[q,qs,score] = self.eval_without_cp_edges(G,x);	
				elseif strcmp(name,'be_continuous')
					[q,qs,score] = self.eval_continuous(G,x);	
				end
			end
		
			if isnan(q);q = 1;end;
			qs = q;
		end
			
		function param = initParam(self,param)
			if ~isfield(param,'numRun');param.numRun = 20;end
			if ~isfield(param,'name') param.name = 'be';end
			if ~isfield(param,'dcor');param.dcor = false;end
			if ~isfield(param,'disp') param.disp = false;end
			
			% for degree correction
			if ~isfield(param,'sigma');param.sigma = 0.1;end
			if ~isfield(param,'beta');param.beta = 0.1;end
			if ~isfield(param,'mu');param.mu = 1e-6;end
			if ~isfield(param,'maxItNum');param.maxItNum = 100;end
			if strcmp( param.name, 'be_cont_dcor' )	
				param.dcor = true;
			end
		end
	end
	
	methods ( Access = private )
		function [C, P, Q, Qs, score] = becont_dcor2( self, G, param )
				
			A = G.adjacency_matrix('binary');
			N = size(A,1);	
			deg = G.degree(); M = G.numEdge();
			
			x0 = rand(N,1);
			x0 = x0/sqrt(x0'*x0);
				
			x = sdpvar( N, 1);
		
			%R = A - deg*deg'/(2*M);	
		
			Obj = -x'*A*x + (deg'*x)^2 / (2*M);
			%Obj = -x'*R*x;
			Const = [x>=0, x'*x <=1];
			
			%options = sdpsettings('verbose',1,'usex0',1);
			options = sdpsettings('verbose',0,'usex0',1);
			%options = sdpsettings('verbose',1);
			assign(x,x0);

			try
				sol = optimize(Const,Obj,options);
				C = value(x);
				P = 1-C;
			catch
				C = ones(N,1)/sqrt(N);
				P = 1-C;	
			end
				
			score = C;
			Q = C'*A*C - (deg'*C)^2 / (2*M);
			Qs = Q;
		end

		function [C, P, Q, Qs, score] = becont_dcor( self, G, param )
			
			A = G.adjacency_matrix('binary');
			N = size(A,1);	
			deg = G.degree(); M = G.numEdge();
			
			x = rand(N,1);x = x/sqrt(x'*x);
			%Qx = (-x'*A*x + ((deg'*x)^2)/(2*M));
			
			%function z = dF(x)
			%	X = x*x'./(1+x*x');
			%	z = sum(X - diag(diag(X)),2) - deg;
			%end
			
			%x = fsolve(@dF, x);	
			%X = x*x'./(1+x*x');
			%R = A - X + diag(diag(X));
			R = (A- deg*deg'/(2*M));	
			
			Qx = -x'*R*x;
			itnum = 0;
			while itnum <= param.maxItNum
				alpha = 1;
					
				while(true)
					%df = param.sigma * (-A* x + (deg)*(deg'*x)/(2*M));
					df = param.sigma * (-R*x/(2*M));
					
					%df = sigma * (Mod * x);
					xnew = x - alpha * df;
					xnew(xnew<0) = 0;
					xnew = xnew/sqrt(xnew'*xnew);
					Qxnew = (-xnew'*R*xnew);
					%Qxnew = (-xnew'*A*xnew + ((deg'*xnew)^2)/(2*M));
					 
					if Qxnew  - Qx  <= param.sigma * df'*(xnew-x) + eps*1e4
						break;
					end
					
					alpha = alpha * param.beta; 
				end
				
				dx = sum(abs(x-xnew));
				if param.disp
					disp(sprintf('itnum=%d, dx=%f',itnum,dx));
				end
				
					
				if dx <= param.mu; 
					x = xnew;
					Qx = Qxnew;
					break;
				end
				x = xnew;
				Qx = Qxnew;
				itnum = itnum + 1;
			end
			C = x;
			P = 1-C;	
			Q = -Qx;
			Qs = Q;
			score = x;	
			
		end
		
		function [Q, Qs, score] = eval_cont_dcor( self, G, x, varargin )
			A = G.adjacency_matrix('binary');
			N = size(A,1);	
			deg = G.degree(); M = G.numEdge();
			Q = (x'*A*x - ((deg'*x)^2)/(2*M));
			Qs = Q;score = x;
		end
		

		function [C, P, Q, Qs, score] = Kernighan_Lin_algorithm_mex( self, G, param )
			
			A = G.adjacency_matrix();
		
			p = nnz(A)/ (size(A,1) * (size(A,1)-1));
			if p > 0.5 % dense network
				A = 1-A;
				A = A - diag(diag(A));
				[r,c,v] = find(triu(A,1));
				P = bealgorithm_without_cp_edges( [r,c], size(A,1), length(r) );	
				C = 1 - P;
				[Q, Qs, score] = self.eval( G, C, P, param );
				
			else % sparse network
				[r,c,v] = find(triu(A,1));
				C = bealgorithm( [r,c], size(A,1), length(r) );	
				P = 1 - C;
				[Q, Qs, score] = self.eval( G, C, P, param );
			end
			
		end
		
		function [C,P,Q,Qs,score] = Kernighan_Lin_algorithm_without_cp_edges_mex(self,G,param)
			A = G.adjacency_matrix('binary');
			[r,c,v] = find( triu(A,1) );
			C = bealgorithm_without_cp_edges( [r,c], size(A,1), length(r) );
			P = 1 - C; 
			[Q, Qs, score] = self.eval( G, C, P, param );
		end
			
		function [C, P, Q, Qs, score] = Kernighan_Lin_algorithm( self, G, param )
			% --------------------------------
			% Initialise
			% --------------------------------
			A = G.adjacency_matrix('binary');
			N = size( A, 1 ); % Number of nodes
			deg = sum( A, 2 ); % deg(i) is the degree of node i
			L = N * ( N - 1  ) / 2; % Number of pairs of nodes 
			M = sum( deg ) / 2; % Number of edges
			while true; C = round( rand( N, 1 ) );if ~all( C ) & ~all( ~C ); break; end;end % Initialise C 
			Ncore = sum( C ); % Number of core nodes
			Mb = N*(N-1)/2 - (N-Ncore)*(N-Ncore-1)/2; % Number of edges in an idealised core-periphery structure 
			Dcore = A * C; % Dcore(i) is the number of edges to core nodes
			[Q, ~, ~] = self.eval( G, C, 1-C, param ); % Calculate the quality of C
			Mcp = C' * deg - C' * Dcore ; % Number of edges between core and peripheral nodes
			x = C;
			
			% --------------------------------
			% Maximise the Borgatti-Everett quality function 
			% --------------------------------
			for j = 1:N		
				fixed = false( N, 1 );  
				Qold = ( Mcp - Mb  *  M / L ) / sqrt( Mb  * ( L - Mb ) );
				dQ = 0;dQmax = - Inf;
				
				for i = 1:ceil( N / 2 )
					% Calculate the increment in the quality obtained by changing node labels 	
					dMb = ( Ncore - N ) .* x + ( N - Ncore - 1 ) .* ( 1 - x ); 
					dMcp = ( Dcore - deg ); dMcp( x==0 ) = -dMcp( x==0 );
					q = ( Mcp + dMcp - ( Mb + dMb )  *  M / L ) ./  sqrt(  ( Mb + dMb ) .* ( L - ( Mb + dMb ) ) );
					q( fixed ) = -Inf; q( isinf( q ) ) = -Inf;
					
					% Update node label to the label that yields the largest increment	
					[qmax, nid] = max( q );
					Ncore = Ncore - ( 2 * x( nid ) - 1 );
					Mb = Mb + dMb( nid );				
					Mcp = Mcp + dMcp( nid );
					Dcore = Dcore + A( :, nid )  *  ( - 2 * x( nid ) + 1 );	
					x( nid ) = 1 - x( nid );
					dQ = dQ + qmax - Qold;
					Qold = qmax;
			
					% Save the core-periphery pair if it attains the largest quality	
					if dQmax < dQ
						xbest = x;
						dQmax = dQ;
						Ncorebest = Ncore; Mbbest = Mb; Mcpbest = Mcp; Dcorebest = Dcore;
					end
					
					fixed( nid ) = true; % Fix the label of node nid
				end
				 	
				if dQmax <= eps;
					break;
				end
				
				Ncore = Ncorebest; Mb = Mbbest; Mcp = Mcpbest; Dcore = Dcorebest;
				x = xbest; C = xbest;
			end
		
			% --------------------------------
			% Calculate the quality of the detected core-periphery pair 
			% --------------------------------
			P = 1 - C; 
			[Q, Qs, score] = self.eval( G, C, P, param );
		end
		
		function [Q, Qs, score] = eval_with_cp( self, G, x, varargin )
			
			if(sum(x)==0)
				Q = 0; Qs =0;score = zeros(G.numNode(),1);
				return;
			end
			
			
			A = G.adjacency_matrix('binary');
			N = size( A, 1 ); 
			
			if nnz(A) == (( N * ( N - 1 ) )) % full matrix
				Q = 1;Qs = 1;score = sparse(N,1);
				return 
			end
			if nnz(A) == 0 % empty matrix
				Q = 0;Qs = 0;score = sparse(N,1);
				return 
			end
		
			Nperi = N-sum( x );
			p = sum( A( : ) )  /  ( N * ( N - 1 ) );
			pb = ( N * ( N - 1 ) -  Nperi * ( Nperi - 1 ) )  /  ( N * ( N - 1 ) );
			vara = p * ( 1 - p ) * N * ( N - 1 ) / 2;	
			varb = pb * ( 1 - pb ) * N * ( N - 1 ) / 2;
			deg = sum( A, 2 ); 
			Dcore = A * x;
			score = ( - deg * pb / 2 + x .* deg / 2 + ( 1 - x ) .* Dcore / 2 ) / ( sqrt( vara ) * sqrt( varb ) );
			score( isnan(score) ) = 0;
			score( isinf(score) ) = 0;
			Q = sum(score);Qs = Q;
			
		end
		
		function [Q, Qs, score] = eval_without_cp_edges( self, G, x, varargin )
			A = G.adjacency_matrix('binary');
			Wcc = x'*A*x;
			Ncore = sum(x);
			N = size(x,1);
			px = Ncore * (Ncore-1)/ (N *(N-1) );
			p = nnz(A) / (N*(N-1));
		
				
			denom = sqrt(p * (1-p) * px *(1-px)) * N*(N-1);
			numer = Wcc - p * Ncore * (Ncore-1);
			Q = numer / denom;
			Qs = Q;score = Q * ones(N,1)/ N;
		end
		
		function [C, P, Q, Qs, score] = becontinuous( self, G, param )
			C = G.max_eigen();
			C = abs(C);P = 1-C;
			[Q,Qs,score] = self.eval( G, C, P, param);
		end
		
		function [Q, Qs, score] = eval_continuous( self, G, x, varargin )
			A = G.adjacency_matrix('binary');
			Q = x'*A*x;Qs = Q;score = x;
		end
	end
end
