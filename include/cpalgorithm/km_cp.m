% Kojaku-Masuda algorithm
%
% @version 0.1.0
% @date 2017/02/19
%
% @author Sadamori Kojaku 
classdef km_cp < cpabst
	
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
			if strcmp( param.null_model, 'erdos-renyi' ) % If the null model is the Erdős–Rényi random graph
				A = G.adjacency_matrix();
				for it = 1:param.numRun 
					ts = cputime;
					if strcmp(param.optalgorithm,'label_switching_matlab')
						[Ct, Pt, Qt, Qst, scoret] = label_switching_potts(self,A,param);
					else
						[Ct, Pt, Qt, Qst, scoret] = self.km_erdos_renyi_mex( G, param ); 
					end
					if Qt > Q
						C = Ct; P = Pt; Q = Qt; Qs = Qst;score = scoret; 
					end
					if param.disp; disp( sprintf( '%d/%d finished ( max Q = %f, %f seconds )', it, param.numRun, full( Q ), cputime - ts ) ); end
				end
			elseif strcmp( param.null_model, 'config' ) % If the null model is the configuratin model 
				A = G.adjacency_matrix();
				[C, P, Q, Qs,score] = self.km_config_mex( G, param ); 
			end
			cpu_time = cputime - cpu_time;
		end
		
		function [Q, Qs, score] = eval( self, G, C, P, param )
			% --------------------------------
			% Calculate Q using the configuration model
			% --------------------------------
			if strcmp( param.null_model, 'config' )
				[Q,Qs,score] = self.eval_km_config(G, C, P, param);
			elseif strcmp(param.null_model,'erdos-renyi')
				[Q,Qs,score] = self.eval_km_erdos_renyi(G, C, P, param);
			else
				Qs = diag(C' * param.Q*C) + 2*diag(C'*param.Q*P);
				Qs =Qs';
				Q = sum(Qs);
				score = Q/size(param.Q,1);
			end
		end
	
		function param = initParam( self, param )
			if isempty( param )	
				param.name = 'km';
			end	
			if ~isfield( param, 'numRun' ); param.numRun = 20; end
			if ~isfield( param, 'optalgorithm' ); param.optalgorithm = 'label_switching'; end
			if ~isfield( param, 'null_model' ); param.null_model = 'erdos-renyi'; end
			if ~isfield( param, 'cptype' ); param.cptype = 'disc'; end
			if ~isfield( param, 'disp' ); param.disp = false; end
			if ~isfield( param, 'K' ); param.K = Inf; end
			if ~isfield( param, 'maxNumCoreNodes' ); param.maxNumCoreNodes = Inf; end
			if ~isfield( param, 'temperature' ); param.temperature = 10; end
			if ~isfield( param, 'isSelfLoop' ); param.isSelfLoop = true; end
		end
	end
	
	methods ( Access = private )
	
		% =========================	
		% KM algorithm (Erdos-renyi)
		% =========================	
		function [C, P, Q, Qs, score] = km_erdos_renyi_mex( self, G, param )
			
			A = G.adjacency_matrix();
			N = size(A,1);
			if isfield(param,'p')
				p = param.p;
			else
				w = sum(A(:));
				if param.isSelfLoop	
					p = w / (N*N);
				else
					p = w / (N*N - N);
				end
			end
			
			
			if param.isSelfLoop	
				[r,c,v] = find(triu(A));
				L = [r,c,v];
			else
				[r,c,v] = find(triu(A,1));
				L = [r,c,v];
			end
		
			if strcmp(param.optalgorithm,'label_switching')
				[cids,x] = km_label_switching( L, size(A,1), size(L,1), min(param.K,N), p, double(param.isSelfLoop) );
			elseif strcmp(param.optalgorithm,'mcmc')
				Nc = min(N,param.maxNumCoreNodes);
				if(Nc==1)
					[cids,x] = starkm_mcmc( L, size(A,1), size(L,1), N*1000, N*100, param.temperature);
				else
					[cids,x] = km_mcmc( L, size(A,1), size(L,1), min(param.K,N), N*1000, N*100, Nc, param.temperature, double(param.isSelfLoop) );
				end
			end
			
			U = sparse(1:N,cids,1);
			C = U;P = U;
			C(x==0,:) = 0;P(x==1,:) = 0;
			
			[Q, Qs, score] = self.eval( G, C, P, param );
		end
		
		function [Q,Qs,score] = eval_km_erdos_renyi( self, G, C, P, param )
			Nc = sum( C ); Np = sum( P ); N = G.numNode();
			if sum(Nc+ Np) ==0
				Q = 0;
				Qs = 0;
				score = 0; 
				return
			end
			A = G.adjacency_matrix();
			if isfield(param,'p')
				p = param.p;
			else
				w = sum(A(:));
				if param.isSelfLoop	
					p = w / (N*N);
				else
					p = w / (N*N - N);
				end
			end
			
			Qs = diag([C+P]'*A*[C+P]) - diag(P'*A*P);
			Qs = Qs - p*( ((Nc + Np).^2)' - (Np.*Np)');
			
			if param.isSelfLoop
				Q = sum(Qs);
				Qs = Qs';
				score = ones(N,1)*Q/N;
			else
				Qs = Qs - ((diag(A) - p)'*C)';
				Qs = Qs';
				Q = sum(Qs);
				score = ones(N,1)*Q/N;
			end
		end
		
		% =========================	
		% Discrete version of core-periphery  structure 
		% (Configuration model)
		% =========================	
		function [C, P, Q, Qs, score] = km_config_mex( self, G, param )	
			A = G.adjacency_matrix();	
			N = size(A,1);
			[r,c,v] = find(triu(A));
			L = [r,c,v];
			
			N = size(A, 1);
			[cids, x, Q, Qs] = km_config_label_switching(L, N, size(L,1), param.numRun);
				
			U = sparse(1:N, cids, 1);
			C = sparse(1:N, cids, x);
			P = U - C;
			score = sparse(N,1,0);	
		end
	
		function [Q,Qs,score] = eval_km_config( self, G, C, P, param )
			if sum(C(:) + P(:)) ==0
				Q = 0;
				Qs = 0;
				score = 0; 
				return
			end
			A = G.adjacency_matrix();
			
			deg = G.degree();
			M = sum(deg)/2;
			U = C + P;	
			Dcore = deg'*C;	
			Dperi = deg'*P;	
			D = deg'*U;
			
			Nc = sum( C ); Np = sum( P ); N = G.numNode();
			if sum(Nc+ Np) ==0
				Q = 0;
				Qs = 0;
				score = 0; 
				return
			end
			
			Qs = diag(U'*A*U) - diag(P'*A*P);
			Qs = Qs - ( (D.^2)' - (Dperi.^2)')/(2*M);
			Qs = Qs/(2*M);
			if param.isSelfLoop
				Q = sum(Qs);
				Qs = Qs';
				score = ones(N,1)*Q/N;
			else
				Qs = Qs - ((diag(A)-deg.^2/(2*M))'*C)';
				Qs = Qs';
				Q = sum(Qs);
				score = ones(N,1)*Q/N;
			end
		end
		
		% =========================	
		% Continuous version of core-periphery  structure 
		% (Configuration model)
		% =========================	
		function [C, P, Q, Qs, score] = km_config_cont_mex( self, G, param )	
			A = G.adjacency_matrix();	
			N = size(A,1);
			[r,c,v] = find(triu(A));
			L = [r,c,v];
			
			if isfield(param,'C') 
				K = size(param.C,2);
				cids = sum((param.C)*sparse(1:K,1:K,1:K),2);
			else
				cids = (1:N)';
			end
			
			[cids,x] = km_config_cont_label_switching( L, N, size(L,1), cids);
			C = diag(x)*sparse(1:N,cids,1);
			P = C;
			[Q, Qs, score] = self.eval( G, C, P, param );
		end
	
		function [Q,Qs,score] = eval_km_config_cont( self, G, C, P, param )
			if sum(C(:)) ==0
				Q = 0;
				Qs = 0;
				score = 0; 
				return
			end
			A = G.adjacency_matrix();
			
			M = sum(A(:))/2;
			d = sum(A,2);
			D = d'*C;	
			W = C'*A*C;
			EW = D'*D/(2*M);
			Qs = diag(W - EW)'/(2*M);
				
			Q = sum(Qs);
			score = ones(size(A,1),1)*Q/size(A,1);
		end

			
		% =========================	
		% Discrete version of core-periphery  structure 
		% (ER random graph model)
		% =========================	
		function [C,P,Q,Qs,score] = label_switching_potts(self,A,param)
			
			% ******************************************* 
			% Main Routine 
			% ******************************************* 
			
			
			% --------------------------------
			% Initialise  
			% --------------------------------
			n = size(A,1);
			%[~,ord] = sort(sum(A,1),'descend');
			active = find(any(A,1));
			
			C = sparse(1:n,1:n,1);P = sparse(n,n,0);
			dQ = 1;
			loopnum = 1;
			iscore = true(n,1);
			Ac = A; Ap = sparse(n,n,0);
			Nc = ones(1,n);Np =zeros(1,n);
			if param.isSelfLoop
				lambda = sum(A(:))/(size(A,1)*size(A,1));	
			else
				lambda = sum(A(:))/(size(A,1)*(size(A,1)-1));	
			end
			while(dQ>0)
				if param.disp;ts = cputime;end
					
				dQ = 0;
				ts = cputime;
				ord = randsample(active,length(active));
				for nid = ord
					% step 1
					cid = find(C(nid,:)+P(nid,:),1,'first');

					cand = find(Ac(nid,:) + Ap(nid,:));
					
					[~,lcid] = find(cid==cand);	
				
					pjgain = (Ac(nid,cand) - lambda* (Nc(cand)-C(nid,cand) ));%vQc;
					if( length(lcid)~=0 & iscore(nid) );pjgain(lcid) = pjgain(lcid) + lambda;end
					cjgain = pjgain + (Ap(nid,cand) - lambda * (Np(cand) - P(nid,cand)));
					if( length(lcid)~=0 & iscore(nid)==false );cjgain(lcid) = cjgain(lcid) + lambda;end
					lgain = 0;
					if length(lcid)~=0
						if iscore(nid)
							lgain = -cjgain(lcid);
						else
							lgain = -pjgain(lcid);
						end
						cjgain(lcid) = -Inf;pjgain(lcid) = -Inf;
					end
					
					% calculate diff of Q vals
					[v2c,v2cid] = max(cjgain);
					[v2p,v2pid] = max(pjgain);
					
					maxhitc = v2c ==cjgain;	
					maxhitp = v2p ==pjgain;	
					if sum(maxhitc)>1 
						maxcand = find(maxhitc);	
						v2cid = randsample(maxcand,1);
					end
					if sum(maxhitp)>1 
						maxcand = find(maxhitp);	
						v2pid = randsample(maxcand,1);
					end
						
					% leave gain
					v2cid = cand(v2cid);
					v2pid = cand(v2pid);
					
					% move 
					%if v2p +lgain <0 & v2c +lgain <0;continue;end
					
					if v2p==v2c & v2c + lgain >0
						if(rand()<0.5)
							v2c = v2c*2;
						else
							v2p = v2p*2;
						end
					end
				
					if v2p < v2c  & v2c + lgain >0% move to core
						if iscore(nid) % move from core
							C(nid,cid) = 0;Nc(cid) = Nc(cid)-1;
							Ac(:,cid) = Ac(:,cid) - A(:,nid);
						else % move from periphery
							P(nid,cid) = 0;Np(cid) = Np(cid)-1;
							Ap(:,cid) = Ap(:,cid) - A(:,nid); 
						end
						
						C(nid,v2cid) = 1;
						Ac(:,v2cid) = Ac(:,v2cid) + A(:,nid); 
						Nc(v2cid) = Nc(v2cid) + 1;
						dQ = dQ + v2c + lgain;
						iscore(nid) = true;
					elseif v2c < v2p & v2p + lgain >0 % move to periphery
						if iscore(nid) % move from core
							C(nid,cid) = 0;Nc(cid) = Nc(cid)-1;
							Ac(:,cid) = Ac(:,cid) - A(:,nid);
						else % move from periphery
							P(nid,cid) = 0;Np(cid) = Np(cid)-1;
							Ap(:,cid) = Ap(:,cid) - A(:,nid); 
						end
						
						P(nid,v2pid) = 1; 
						Ap(:,v2pid) = Ap(:,v2pid) + A(:,nid); 
						Np(v2pid) = Np(v2pid) + 1;
						dQ = dQ + v2p + lgain;
						iscore(nid) = false;
					end
				end
				
				
				remove = Nc ==0 & Np ==0;
				if any(remove) 
					C(:,remove)=[];P(:,remove)=[];
					Nc(remove)=[];Np(remove)=[];
					Ac(:,remove)=[];Ap(:,remove)=[];
				end

				%if param.disp;dt = cputime-ts;disp(sprintf('   %dth loop took %f seconds dQ = %f',loopnum,dt,full(dQ)));end
				loopnum = loopnum + 1;
				dt =cputime-ts;
			end
			[Q,Qs,score] = self.eval_km_erdos_renyi(graph(A), C, P, param);
		end	
		
		function [C,P,Q,Qs,score] = label_switching(self,M,model)
			
			% ******************************************* 
			% Main Routine 
			% ******************************************* 
			
			% --------------------------------
			% Initialise  
			% --------------------------------
			n = size(M,1);
			P = sparse(1:n,1:n,1);C = sparse(n,n,0);
			
			% --------------------------------
			% Label switching (Main routine) 
			% --------------------------------
			n = size(M,1);
			dQ = 1;loopnum = 1;
			while(dQ>1e-6)
				dQ = 0;
				ord = randsample(n,n)';
				for nid =ord
					iscore = any(C(nid,:));
					cid = find(C(nid,:) + P(nid,:)>0);
					cQv = (C'*M(:,nid))';	
					pQv = (P'*M(:,nid))';	
					vQc = M(nid,:)*C;
					vQp = M(nid,:)*P;
					vQv = M(nid,nid);
				
					if iscore
						cQv(cid) = cQv(cid)-vQv;
						vQc(cid) = vQc(cid)-vQv;
					else
						pQv(cid) = pQv(cid)-vQv;
						vQp(cid) = vQp(cid)-vQv;
					end 
					
					% core join gain 
					cjgain = cQv  + vQc+ pQv + vQp + vQv;
					
					% periphery join gain 
					pjgain = cQv + vQc ;
					% calculate diff of Q vals
					[v2c,v2cid] = max(cjgain);
					[v2p,v2pid] = max(pjgain);
					
					if iscore
						% leave gain
						cid = find(C(nid,:)>0);
						lgain = -cQv - pQv - vQc - vQp - vQv;
						lgain = lgain(cid);
					else
						% leave gain
						cid = find(P(nid,:)>0);
						lgain = -cQv - vQc;
						%lgain = -cQv - pQv - vQc - vQp + vQp + pQv;
						lgain = lgain(cid);
					end
					
				
					if v2p < v2c  & v2c + lgain >0% move to core
						if iscore % move from core
							C(nid,cid) = 0;
						else % move from periphery
							P(nid,cid) = 0;
						end
						
						C(nid,v2cid) = 1;
						dQ = dQ + v2c + lgain;
					elseif v2c < v2p & v2p + lgain >0 % move to periphery
						if iscore % move from core
							C(nid,cid) = 0;
						else % move from periphery
							P(nid,cid) = 0;
						end
						
						P(nid,v2pid) = 1; 
						dQ = dQ + v2p + lgain;
					end
				end
				remove = ~any(C+P,1);
				if any(remove) 
					C(:,remove)=[];P(:,remove)=[];
				end
				loopnum = loopnum + 1;
			end
			
			slice = any(C+P,1);
			C = C(:,slice);
			P = P(:,slice);
			Nc = sum(C,1); Np = sum(P,1);
			score = sum( (M*(C+P)).*C,2);
			Qs = diag( (C+P)'*M*(C+P) - P'*M*P)';
			[Qs,ord] = sort(Qs,'descend');
			C = C(:,ord);P = P(:,ord);
			Q = sum(Qs);
		end
	end
end
