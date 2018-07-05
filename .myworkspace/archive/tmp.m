
		function [C,P,Q,Qs] = label_switching_potts(self,A,model)
			
			% ******************************************* 
			% Main Routine 
			% ******************************************* 
			
			% --------------------------------
			% Initialise  
			% --------------------------------
			
			% --------------------------------
			% Label switching (Main routine) 
			% --------------------------------
			n = size(A,1);
			C = sparse(1:n,1:n,1);P = sparse(n,n,0);
			dQ = 0; 
				
			active = find(any(A,1));
			ord = randsample(active,length(active));
			
			dQ = 1;
			loopnum = 1;
			Nc = ones(1,n);Np =zeros(1,n);
			Ac = A; 
			Ap = sparse(n,n,0);
			iscore = true(n,1);
			while(dQ>0 & loopnum <= model.maxitnum)
				dQ = 0;
				if model.disp
					ts = cputime;
				end
				ts =cputime;
				
				% need update
				for nid = ord
					% step 1
					cid = find(C(nid,:)+P(nid,:),1,'first');	
					%vQc = A(nid,:)*C - model.lambda * (Nc-C(nid,:));
					%vQp = A(nid,:)*P - model.lambda * (Np-P(nid,:));
					
					cand = find(Ac(nid,:)+Ap(nid,:));
						
					lid = cand==cid;
					vQc = Ac(nid,cand) - model.lambda * Nc(cand);vQc(lid) = vQc(lid) - model.lambda;
					vQp = Ap(nid,cand) - model.lambda * Np(cand);vQp(lid) = vQp(lid) - model.lambda;
					
					% periphery join gain 
					pjgain = 2*vQc;
					
					% core join gain 
					cjgain = pjgain + 2*vQp + A(nid,nid);
				
					
					% calculate diff of Q vals
					[v2c,v2cid] = max(cjgain);
					[v2p,v2pid] = max(pjgain);
					
					if cid==v2cid;continue;end	
					v2cid = cand(v2cid);
					v2pid = cand(v2pid);
					
					if iscore(nid)
						% leave gain
						lgain = -cjgain(lid) + 2*A(nid,nid);
					else
						% leave gain
						lgain = -2*vQc(lid);
					end
					
					if v2p < v2c & 0 <= v2c + lgain % move to core
						
						if iscore(nid)
							C(nid,cid) = 0;Nc(cid) = Nc(cid)-1;
						else
							P(nid,cid) = 0;Np(cid) = Np(cid)-1;
						end
						C(nid,v2cid) = 1;
						cids(nid) = v2cid;
						 
						Ac(:,cid) = Ac(:,cid) - A(:,nid); 
						Ac(:,v2cid) = Ac(:,v2cid) + A(:,nid); 
							
						Nc(v2cid) = Nc(v2cid) + 1;
						dQ = dQ + v2c + lgain;
						
					elseif v2c <= v2p & 0 <= v2p + lgain % move to periphery
						if iscore(nid)
							C(nid,cid) = 0;Nc(cid) = Nc(cid)-1;
						else
							P(nid,cid) = 0;Np(cid) = Np(cid)-1;
						end
						P(nid,v2pid) = 1; 
						cids(v2pid) = 1;

						Ap(:,cid) = Ap(:,cid) - A(:,nid); 
						Ap(:,v2cid) = Ap(:,v2cid) + A(:,nid); 

						Np(v2pid) = Np(v2pid) + 1;
						dQ = dQ + v2p + lgain;
					end
				end
				remove = Nc ==0 & Np ==0;
				if any(remove)
					C(:,remove)=[];P(:,remove)=[];
					Nc(remove)=[];Np(remove)=[];
					Ac(:,remove)=[];Ap(:,remove)=[];
				end
				dt =cputime-ts;
				dt

				if model.disp
					dt = cputime-ts;
					disp(sprintf('   %dth loop took %f seconds dQ = %f',loopnum,dt,dQ))	
				end
				loopnum = loopnum + 1;
			end
			slice = sum(C,1)==0 ;
			C(:,slice)=[];P(:,slice)=[];
			Nc(slice)=[];Np(slice)=[];
			
			Nc = sum(C,1); Np = sum(P,1);
			Qs = diag( (C+P)'*A*(C+P) - P'*A*P)' - model.lambda * (Nc+Np).^2  + model.lambda * Np.^2;
			Qs = Qs / sum(A(:));
			Q = sum(Qs);
			
			
		end	
