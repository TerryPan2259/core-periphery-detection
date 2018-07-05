% Discrete variant of the Borgatti-Everett algorithm 
%
% @version 0.1.0
% @date 2017/02/19
%
% @author Sadamori Kojaku 
classdef cucuringu_cp < cpabst

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
			if strcmp( param.name, 'lap_core')
				[C,P,Q,Qs,score] = self.lap_core(G, param);
			elseif strcmp( param.name, 'lapsgn_core')
				[C,P,Q,Qs,score] = self.lapsgn_core(G, param);
			elseif strcmp( param.name, 'low_rank_core')
				[C,P,Q,Qs,score] = self.low_rank_core(G, param);
			end
			
			cpu_time = cputime-ts;
		end
		
		function [q,qs,score] = eval(self,G,x,varargin)
		
			N = G.numNode();
			A = G.adjacency_matrix();
			Nc = sum(x);	
			
			Mcc = x'*A*x/2; 
			Mcp = x'*A*(1-x);
			Mpp = (1-x)'*A*(1-x);
			
			q = Mcc/(Nc*(Nc-1)/2) + Mcp/(Nc*(N-Nc)) - Mpp / ((N-Nc)*((N-Nc)-1)/2);
		
			if nargin ==5
				param = varargin{2};
				q = q - param.gamma * (  Nc/N - param.beta   );
			end
			qs = q;
			score = x;
		end
			
		function param = initParam(self,param)
			if ~isfield(param,'name') param.name = 'low_rank_core';end
			if ~isfield(param,'disp') param.disp = false;end
			if ~isfield(param,'beta') param.beta = 0.1;end
			if ~isfield(param,'gamma') param.gamma = 0;end
		end
	end
	
	methods ( Access = private )	
		
	
		function [C,P,Q] = find_cut(self, G, score, b)
			[val,ord] = sort(score,'descend');
			N = G.numNode();
			A = G.adjacency_matrix();
			q = zeros(N,1);
			deg = G.degree();M = G.numEdge();N = G.numNode();
			for i = b:N-b
				Mcc = sum(sum( A( ord(1:i),ord(1:i) ) ))/2;
				Mcp = sum(sum( A( ord(1:i),ord(i+1:end) ) ));
				Mpp = sum(sum( A( ord(i+1:end),ord(i+1:end) ) ))/2;
				q(i) = Mcc/(i*(i-1)/2) + Mcp/(i*(N-i)) - Mpp / ((N-i)*((N-i)-1)/2);
			end
			
			[Q,idx] = max(q);
			Q = Q/N;
			C  = sparse( ord(1:idx),1,1,N,1);
			P  = 1-C;
		end
		
		function [C,P,Q,Qs,score] = lapsgn_core(self, G, param)
			T = G.transition_matrix();
			N = size(T,1);
			try
				[v,d] = eigs(T+eye(N), 1,'sm');
			catch
				[v,d] = eig(full(T));
				[val,ord] = sort(diag(d));
				v = v(:, ord(1));
				d = d(ord(1));
			end
			v = sign(v);
			C = sparse(N,1,0);C(v>0) = 1;
			Q = self.eval( G, C, 1-C, param); 
			q1 = self.eval( G, 1-C, C, param);
			score = v;	
			if q1 > Q
				score = -score;
				C = 1-C;Q = q1; 
			end
			Qs = Q;
			P = 1-C; 
		end
			
		
		function [C,P,Q,Qs,score] = lap_core(self, G, param)
			T = G.transition_matrix();
			N = size(T,1);
			try
				[v,d] = eigs(T+eye(N), 1,'sm');
			catch
				[v,d] = eig(full(T));
				[val,ord] = sort(diag(d));
				v = v(:, ord(1));
				d = d(ord(1));
			end
			
			score = v;
			[C,P,Q] = self.find_cut(G,score,ceil(param.beta*N));	
			Qs = Q;
			
		end
		
		function [C,P,Q,Qs,score] = low_rank_core(self, G, param)
			
			A = G.adjacency_matrix();
			try
				[v,d] = eigs(A,2,'lm');
			catch
				[v,d] = eig(full(A));
				[val,ord] = sort(diag(d),'descend');
				v = v(:, ord(1:2));
				d = diag(d(ord(1:2)));
			end	
			N = size(A,1);
			At = v*d*v';
			At(At>0.5) = 1;	
			At(At<=0.5) = 0;
			score = sum(At,2);
			[C,P,Q] = self.find_cut(G,score,ceil(param.beta*N));	
			Qs = Q;
		end
			
	end
end
