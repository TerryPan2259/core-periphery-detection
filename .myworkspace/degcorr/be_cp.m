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
						[Ct,Pt,Qt,Qst,scoret] = self.becont_dcor(G,param);
					end
							
					if Qt > Q | it == 1
						C = Ct; P = Pt; Q = Qt; score = scoret; Qs = Qst;
					end
				end
			else
				for it = 1:param.numRun
					if strcmp(param.name,'be') 
						[Ct,Pt,Qt,Qst,scoret] = self.Kernighan_Lin_algorithm_mex(G,param);
					elseif strcmp(param.name,'be_diag')
						[Ct,Pt,Qt,Qst,scoret] = self.Kernighan_Lin_algorithm_without_cp_mex(G,param);
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
				elseif strcmp(name,'be_diag')
					[q,qs,score] = self.eval_without_cp(G,x);	
				elseif strcmp(name,'be_continuous')
					[q,qs,score] = self.eval_continuous(G,x);	
				end
			end
		end
			
		function param = initParam(self,param)
			if ~isfield(param,'numRun');param.numRun = 20;end
			if ~isfield(param,'name') param.name = 'be';end
			if ~isfield(param,'dcor');param.dcor = false;end
			if ~isfield(param,'disp') param.disp = false;end
			
			% for degree correction
			if ~isfield(param,'sigma');param.sigma = 0.1;end
			if ~isfield(param,'beta');param.beta = 0.1;end
			if ~isfield(param,'mu');param.mu = 1e-4;end
			if strcmp( param.name, 'be_cont_dcor' )	
				param.dcor = true;
			end
		end
	end
	
	methods ( Access = private )

		function [C, P, Q, Qs, score] = becont_dcor( self, G, param )
			
			A = G.adjacency_matrix('binary');
			N = size(A,1);	
			deg = G.degree(); M = G.numEdge();
			
			x = rand(N,1);x = x/sqrt(x'*x);
			Qx = (-x'*A*x + ((deg'*x)^2)/(2*M));
			while(true)
				alpha = 1;	
				while(true)
					df = param.sigma * (-A* x + (deg)*(deg'*x)/(2*M));
					%df = sigma * (Mod * x);
					xnew = x - alpha * df;
					xnew(xnew<0) = 0;
					xnew = xnew/sqrt(xnew'*xnew);
					Qxnew = (-xnew'*A*xnew + ((deg'*xnew)^2)/(2*M));
					if Qxnew  - Qx  < param.sigma * df'*(xnew-x)
						break;
					end
					alpha = alpha * param.beta; 
				end
				if abs(Qxnew  - Qx)/N < param.mu; 
					x = xnew;
					Qx = Qxnew;
					break;
				end
				x = xnew;
				Qx = Qxnew;
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
			x = rand(N,1);x = x/sqrt(x'*x);
			Q = (-x'*A*x + ((deg'*x)^2)/(2*M));
			Qs = Q;score = x;
		end
		

		function [C, P, Q, Qs, score] = Kernighan_Lin_algorithm_mex( self, G, param )
			A = G.adjacency_matrix();	
			[r,c,v] = find(triu(A,1));
			C = bealgorithm( [r,c], size(A,1), length(r) );	
			P = 1 - C;
			[Q, Qs, score] = self.eval( G, C, P, param );
		end
		
		function [C,P,Q,Qs,score] = Kernighan_Lin_algorithm_without_cp_mex(self,G,param)
			A = G.adjacency_matrix('binary');
			[r,c,v] = find( triu(A,1) );
			C = bealgorithm_diag([r,c],G.numNode(),length(r));
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
		
		function [C,P,Q,Qs,score] = Kernighan_Lin_algorithm_without_cp(self,G,param)
			A = G.adjacency_matrix('binary');
			n = size(A,1);
			deg = sum(A,2);
			M = n*(n-1)/2;
			muy = sum(deg)/2;
			
			Q = 0;
			% init x
			if n<=2
				C = ones(n,1);P=1-C;Q=1;score = zeros(n,1);
				return;
			end

			while(true)
				x = round(rand(n,1));
				if ~all(x) & ~all(~x);break;end
			end
		
			C = x;
			% update x
			for j = 1:n
				assigned =false(n,1);
				nc = sum(x);np = n-nc;
				cedges = x'*A*x/2;
				pedges = (1-x)'*A*(1-x)/2;
				Ec = cedges;	
				Ep = pedges;	
				pa = (Ec + Ep)/(nc*(nc-1)/2  + np*(np-1)/2);pb = (nc*(nc-1)/2)/(nc*(nc-1)/2  + np*(np-1)/2);
				Qold = (1-pa)*Ec - ((nc-1)*nc/2-Ec)*pa;
				Qold = Qold/( sqrt(pa*(1-pa)*pb*(1-pb)) );
				indeg = A*x;
				if isinf(Qold);Qold = 0;end;
				xnum = sum(x);	
	
				dQ = 0;dQmax = -Inf;
				for i = 1:ceil(n/2)
					% compute score
					score = sparse(n,1,0);
					nc = xnum +(1-2*x);np = n-nc;
					Ec = cedges + (1-2*x).*indeg;
					Ep = pedges + (1-2*(1-x)).*(deg-indeg);
					pa = (Ep+Ec)./(nc.*(nc-1)/2  + np.*(np-1)/2);
					pb = (nc.*(nc-1)/2)./(nc.*(nc-1)/2  + np.*(np-1)/2);
					score = (1-pa).*Ec - ((nc-1).*nc/2-Ec).*pa;
					score = score./( sqrt(pa.*(1-pa).*pb.*(1-pb)) );
					score(isinf(score)) = 0;
					score(isnan(score)) = 0;
				
					
					[qmax,idx] = max(score);
					assigned(idx) = true;
					cedges = cedges + (1-2*x(idx))*indeg(idx);
					pedges = pedges + (1-2*(1-x(idx)))*(deg(idx)-indeg(idx));
					indeg = indeg + (1-2*x(idx))*A(:,idx);
					xnum = xnum + (1-2*x(idx));
					x(idx) = 1-x(idx);
					dQ = dQ + qmax - Qold;
					Qold = qmax;
						
					if dQmax < dQ
						xbest = x;
						dQmax = dQ;
					end
				end
					
				if dQmax <=eps*10000;
					break;
				end
				x = xbest;
				C = xbest;
			end
			P = 1-C;[Q,Qs,score] = self.eval(G, C, P, param);	
		end
		
		function [Q, Qs, score] = eval_with_cp( self, G, x, varargin )
			if(sum(x)==0)
				Q = 0; Qs =0;score = zeros(G.numNode(),1);
				return;
			end
			A = G.adjacency_matrix('binary');
			N = size( A, 1 ); 
			Nperi = N-sum( x );
			p = sum( A( : ) )  /  ( N * ( N - 1 ) );
			pb = ( N * ( N - 1 ) -  Nperi * ( Nperi - 1 ) )  /  ( N * ( N - 1 ) );
			vara = p * ( 1 - p ) * N * ( N - 1 ) / 2;	
			varb = pb * ( 1 - pb ) * N * ( N - 1 ) / 2;
			deg = sum( A, 2 ); 
			Dcore = A * x;
			score = ( - deg * pb / 2 + x .* deg / 2 + ( 1 - x ) .* Dcore / 2 ) / ( sqrt( vara ) * sqrt( varb ) );
			Q = sum(score);Qs = Q;
		end
		
		function [Q, Qs, score] = eval_without_cp( self, G, x, varargin )
			A = G.adjacency_matrix('binary');
			n = size(A,1);
			cores = find(x);	
			peripheries = find(x==0);	
			L = [];
			if length(cores)>1
				t = nchoosek(cores,2);
				L =[L;t];
			end
			if length(peripheries)>1
				t = nchoosek(peripheries,2);
				L =[L;t];
			end
			
			eids = sub2ind([n,n],L(:,1),L(:,2));
			a = A(eids);	
			xx = x(L(:,1)) + x(L(:,2)) - ( x(L(:,1)) .* x(L(:,2)) );	
			Q = corr(a,xx);	
			
			if isnan(q)
				Q = 0;
			end
			Qs = Q;
			score = sum(A,2);	
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
