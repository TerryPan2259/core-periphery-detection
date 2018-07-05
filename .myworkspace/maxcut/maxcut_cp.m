% @version 0.1.0
% @date 2017/02/19
%
% @author Sadamori Kojaku 
classdef maxcut_cp < cpabst
	
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
			P =[];C = []; Q = - Inf;Qs = [];score = [];
			cpu_time = cputime;
	
			% --------------------------------
			% Optimise node labels using a label switching algorithm 
			% --------------------------------
			for it = 1:param.numRun 
				ts = cputime;
				if strcmp(param.optalgorithm, 'mcmc')
					[Ct, Pt, Qt, Qst,scoret] = self.mcmc_mex( G, param ); 
				elseif strcmp(param.optalgorithm, 'kl')
					[Ct, Pt, Qt, Qst,scoret] = self.kl_mex( G, param ); 
				end
				if Qt > Q
					C = Ct; Q = Qt; P = Pt;Qs = Qst;score = scoret;
				end
				if param.disp; disp( sprintf( '%d/%d finished ( max Q = %f, %f seconds )', it, param.numRun, full( Q ), cputime - ts ) ); end
			end
			cpu_time = cputime - cpu_time;
		end
		
		function [Q,Qs,score] = eval( self, G, C, P, param )
			A = G.adjacency_matrix();
			N = size(A,1);
			for i = 1:size(C,2)
				Qs(1,i) = C(:,i)'*A*P(:,i);
			end
			Q = sum(Qs);
			score = Q*ones(N,1)/N;
		end
	
		function param = initParam( self, param )
			if isempty( param )	
				param.name = 'maxcut';
			end	
			if ~isfield( param, 'numRun' ); param.numRun = 20; end
			if ~isfield( param, 'optalgorithm' ); param.optalgorithm = 'mcmc'; end
			if ~isfield( param, 'temperature' ); param.temperature = 10; end
			if ~isfield( param, 'disp' ); param.disp = false; end
		end
	end
	
	methods ( Access = private )
	
		% =========================	
		% MCMC algorithm 
		% =========================	
		function [C, P, Q, Qs,score] = mcmc_mex( self, G, param )
				
			A = G.adjacency_matrix();	
			N = size(A,1);
			[r,c,v] = find(triu(A,1));
			L = [r,c,v];
			[S,Q] = maxcut_mcmc( L, size(L,1), double(rand(N,1)<0.5), N, N*1000, N*100, param.temperature );
			K = max(S(:,3));
			C = sign(sparse( S(:,1), S(:,3),1,N,K));
			P = sign(sparse( S(:,2), S(:,3),1,N,K));
			[Q,Qs,score] = self.eval( G, C, P, param );
		end
		
		% =========================	
		% Kernighan-Lin algorithm 
		% =========================	
		function [C, P, Q, Qs,score] = kl_mex( self, G, param )
				
			A = G.adjacency_matrix();	
			N = size(A,1);
			[r,c,v] = find(triu(A,1));
			L = [r,c,v];
			[S,Q] = maxcut_kl( L, size(L,1), double(rand(N,1)<0.5), N);
			K = max(S(:,3));
			C = sign(sparse( S(:,1), S(:,3),1,N,K));
			P = sign(sparse( S(:,2), S(:,3),1,N,K));
			[Q,Qs,score] = self.eval( G, C, P, param );
		end
	end
end
