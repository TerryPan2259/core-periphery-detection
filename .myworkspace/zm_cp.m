% Zhang et.al core-periphery detection algorithm 
%
% @version 0.1.0
% @date 2017/02/19
%
% @author Sadamori Kojaku 
classdef zm_cp < cpabst

	methods ( Access = public )

		function [C, P, Q, Qs, score, param, cpu_time] = detect( self, G, param )
			% --------------------------------
			% Initialise
			% --------------------------------
			C= []; P = []; Q = -Inf; score = [];
			
			% --------------------------------
			% Optimise node labels using the Kernighan-Lin algorithm 
			% --------------------------------
			ts = cputime;
			for it = 1:param.numRun
				if strcmp( param.optalgorithm, 'em-algorithm' ) % If the null model is the Erdős–Rényi random graph
					[Ct,Pt,Qt,Qst,scoret] = self.EMalgorithm_mex(G,param);
				elseif strcmp( param.optalgorithm, 'spectral-clustering' )	
					[Ct,Pt,Qt,Qst,scoret] = self.spectral_clustering(G,param);
				end
				if Qt > Q
					C = Ct; P = Pt; Q = Qt; score = scoret; Qs = Qst;
				end
			end
			cpu_time = cputime-ts;
		end
		
		function [Q,Qs,score] = eval(self,G,x,varargin)
			A = G.adjacency_matrix('binary');
			U = [x,1-x];
			W = U'*A*U;
			deg = sum(A,2);
			Ddeg = (deg'*U)';
			th = sum( diag(deg) * U * diag(1./Ddeg),2);
			
			score = gammaln(deg + 1);
			score(isinf(score)) = 0;
			score(isnan(score)) = 0;
			L = W.*log( W./(Ddeg*Ddeg') );
			L(isnan(L))=0;
			L(isinf(L)) = 0;
			Q = sum(deg)/2 + sum(score)/2 + sum(L(:))/2;			
			Qs = Q;
		end
			
		function param = initParam(self,param)
			if ~isfield(param,'numRun');param.numRun = 20;end
			if ~isfield(param,'name') param.name = 'zm';end
			if ~isfield(param,'disp') param.disp = false;end
			if ~isfield(param,'maxItNum_bp') param.maxItNum_bp = 100;end
			if ~isfield(param,'maxItNum_em') param.maxItNum_em = 100;end
			if ~isfield(param,'tol') param.tol = 1e-4;end
			if ~isfield(param,'optalgorithm') param.optalgorithm = 'em-algorithm';end
		end
	end
	
	methods ( Access = private )
		
			
		function [C,P,Q,Qs,score] = spectral_clustering( self, G, param )
			A = G.adjacency_matrix('binary');
			[v,d] = eigs(A,1,'lm');
			X = v * d;
			labs = kmeans(X,2);
			C = sparse(1:length(labs),labs,1);
			
			E = C'*A*C;
			Ns = sum(C,1);
			
			if E(1,1)/(Ns(1)*(Ns(1)-1)) < E(2,2)/(Ns(2)*(Ns(2)-1));
				C = C(:,2);
			else
				C = C(:,1);
			end
			[Q,Qs,score] = self.eval( G, C);
			P = 1-C;	
		end
		
		function [C,P,Q,Qs,score] = EMalgorithm_mex( self, G, param )
			A = G.adjacency_matrix('binary');
			[r,c,v] = find(triu(A));
			while true
				C = zm_bp( [r,c,v], size(A,1), length(r));
				C = double(C==1);
				P = 1-C;
			
				if any(C)== false | all(C)
				else
					break;
				end
			end
	
			D = sum([C,P]' * A,2);
			N = sum([C,P],1);
			if D(1)/N(1) < D(2)/N(2)
				C = 1-C;
				P = 1-P;
			end
			[Q,Qs,score] = self.eval( G, C);
		end
		function [C,P,Q,Qs,score] = EMalgorithm( self, G, param )
		
			% ============================	
			% Initialise
			% ============================	
			N = G.numNode();
			A = G.adjacency_matrix('binary');
			[rnodes,cnodes] = find(A);
			prr = rand();prs = rand();pss = rand();
			p = [prr,prs;prs,pss]; 	
			gamma = rand(2,1);gamma = gamma / sum(gamma);
			q = rand(N,2);
				
			% ============================	
			% EMalgorithm 
			% ============================
			while(true)
				oldp = p;
				oldgamma = gamma;
				
				% E-step	
				[q,eta] = self.brief_propagation( p, gamma, q, A, param );
				
				% M-step
				gamma = sum(q,1)/N;
				qrsij = zeros( length(rnodes),4);
				rs = [1,1;1,2;2,1;2,2];	
				for mid = 1:length(rnodes)
					i = rnodes(mid);
					j = cnodes(mid);
					for l = 1:size(rs,1)
						r = rs(l,1);s = rs(l,2);
						idx = find( rnodes == j & cnodes == i); 
						qrsij(mid,l) = qrsij(mid,l) + eta(mid, r) * eta(idx,s) * p(r,s);
					end
					qrsij(mid,:) = qrsij(mid,:) / sum(qrsij(mid,:)); 
				end
				numer = sum( qrsij );
				for l = 1:size(rs,1)
					r = rs(l,1);s = rs(l,2);
					p(r,s) = numer(l) / ( sum(q(:,r) ) * sum(q(:,s)));		
				end
				
				dif = sum(abs( gamma(:) - oldgamma(:)))/2 + sum(abs( p(:) - oldp(:)))/4;
				if dif < param.tol
					break
				end	
			end
			
			[val,cid] = max(q');
			U = sparse( 1:N, cid,1,N,2);
			if p(2,2) < p(1,1)
				C = U(:,1); P = U(:,2);
			else
				C = U(:,2); P = U(:,1);
			end
			[Q,Qs,score] = self.eval( G, C);
		end	
		
		function [q, eta] = brief_propagation( self, p, gamma, q, A, param )
			N = size(A,1);
			[r,c] = find(A);
			L = [r,c];
			
			eta = rand( size(L,1), 2 ); % initialise
			
			% update eta
			itnum = 0;
			while( itnum < param.maxItNum )
				itnum = itnum +1;
				neweta = eta;
				for mid = 1:size(L,1)
					i = L(mid,1);
					j = L(mid,1);
					for r = 1:2
						pd = 1;
						for k = 1:size(A,2)
							if A(i,k)~=0;continue;end
							pd = pd * ( 1  - q(k,:)*p( r,:)'  );
						end
						
						for k = 1:size(A,2)
							if A(i,k)==0;continue;end
							if k==j;continue;end
							idx = find(((L(:,1) == k) & (L(:,2) == i)));
							pd = pd * ( neweta( idx,:)*p(r,:)'  );
							%pd = pd * ( eta( idx,:)*p(r,:)'  );
						end
						neweta(mid,r) = pd * gamma(r);
					end
					neweta(mid,:) = neweta(mid,:) / sum(neweta(mid,:));
				end
				% update q
				newq = q;
				for i = 1:N
					for r = 1:2
						pd = 1;
						for k = 1:size(A,2)
							if A(i,k)~=0;continue;end
							pd = pd * ( 1  - newq(k,:)*p(r,:)'  );
							%pd = pd * ( 1  - q(k,:)*p(r,:)'  );
						end
						for k = 1:size(A,2)
							if A(i,k)==0;continue;end
							idx = find(((L(:,1) == k) & (L(:,2) == i)));
							pd = pd * ( neweta( idx,:)*p(r,:)' );
						end
						newq(i,r) = gamma(r) * pd;	
					end	
					newq(i,:) = newq(i,:) / sum(newq(i,:));
				end
				
				dif = sum(sum(abs(neweta - eta)))/( size(L,1)*2 ) + sum(sum((abs(newq-q))))/( N*2 );
				
				if dif < param.tol 
					break;	
				end	
				eta = neweta;
				q = newq;
			end
		end	
		
	end
end
