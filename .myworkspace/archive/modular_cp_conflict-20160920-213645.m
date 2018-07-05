classdef modular_cp 
	methods (Access = public)
		function [C,P,Q,score] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
		
			switch model.solver
				case 'louvain'
					[C,P,Q,score] = self.louvain(graph(A),model);
				otherwise	
					disp('unknown solver for borgatt & Everett model')
			end
		end
	end
	methods (Access = private)
		function [C,P,Qs,score] = louvain(self,G,model)
			
			C = []; P = [];Qs = [];score = []; Q = -Inf;
			for it = 1:model.itnum	
				[Ct,Pt,Qt,Qst,scoret] = self.louvain_one_run(G,model);
				
				if Qt > Q
					C = Ct;P = Pt;Q = Qt;Qs = Qst;score = scoret; 
				end
			end
		end
			
		function [C,P,Q,Qs,score] = louvain_one_run(self,G,model)
			
			% ******************************************* 
			% Main Routine 
			% ******************************************* 
			
			% --------------------------------
			% Initialise  
			% --------------------------------
			T = G.transition_matrix();
			d = G.degree();	
			n = G.sz();
			m = G.edgenum();
			st = d/(2*m);
			
			% --------------------------------
			% Node Penality  
			% --------------------------------
			%q = model.lambda * (sum(T,2)/n - d/(2*m));
			q = model.lambda * (sum(T,2) - ones(n,1))/n;
				
			% --------------------------------
			% Louvain Algorithm (Main routine) 
			% --------------------------------
			C = sparse(1:n,1:n,1);
			dQ = 0; qt = q;
			Gt = G;
			while(true)
				[Ct,dq] = first_step(Gt,qt);	
				if dq==0;break;end
				dQ = dQ + dq;
					
				qt = Ct'*qt;
				C = C *Ct;	
				Gt = second_step(Gt,Ct);
			end
			if size(C,2)==0
				C = sparse(G.sz(),1,0);
				P = 1-C;
				score = C;
				Qs = [0];Q =0;
				return
			end
			
			score = calc_node_score(C,G,q);
			Qs = score'*C;Q = sum(Qs);
			[Qs,ord] = sort(Qs,'descend');
			C = C(:,ord);	
			T = G.transition_matrix();
			P = find_periphery(G,C);
			
			P(:,sum(C,1)<=1)=[];
			C(:,sum(C,1)<=1)=[];

	
			% ******************************************* 
			% Subroutine 
			% ******************************************* 
			
			function [U,dQ] = first_step(G,q)
				dQ = 0;
				d = G.degree();
				A = G.adjacency_matrix();
				T = G.transition_matrix();
				n = G.sz();
				m = G.edgenum();	
				U = sparse(1:n,1:n,1);
				dQ = 0;
				vAv = diag(A);	
				
				M = (A - d*d'/(2*m))/(2*m); 
				%M = (T  - ones(n)/n)/n;
				%M = (M + M')/2;
				
				vMv = diag(M);
				M = M - diag(vMv);
 
				[~,ord] = sort(rand(1,n));
				update = true;
				while(update)
					update = false;
					for nid = ord
						iscore = any(U(nid,:),2);
						isbelonged = any(U(nid,:),1);
						if iscore
							vAc = A(nid,:)*U - vAv(nid)*U(nid,:);
							dU = d'*U;dU(isbelonged) = dU(isbelonged) - d(nid);
							
							qp = M(nid,:)*U;
							qr = -qp(isbelonged) -  q(nid) - vMv(nid);
							qp = qp - qp(isbelonged);
						
							%qr = -2*q(nid) - ( (vAc(isbelonged) + vAv(nid) - d(nid)*(d(nid)+dU(isbelonged))/(2*m))/(2*m)  ); 
							%qp = vAc/(2*m) -d(nid)*dU/(4*m^2);
						else
							vAc = A(nid,:)*U;
							dU = d'*U;
							qr = 0;
								
							qp = M(nid,:)*U + q(nid) + vMv(nid);	
							%qp = vAc/(2*m) -d(nid)*dU/(4*m^2) + 2*q(nid);
							%qp = qp + (vAv(nid)/(2*m) - d(nid)^2/(4*m^2));
						end
						%qp = vAc/(2*m) -d(nid)*dU/(4*m^2) + 2*q(nid)*ones(1,size(u,2));
						%qp = (qp + qr * ones(1,size(U,2))).*~any(U(nid,:),1);
						[dq,cid] = max([qp,qr]);
						if dq <=eps*10000;continue;end
						if cid <=size(U,2)
							if isbelonged(cid)==1;continue;end
							update = true;
							dQ = dQ + 2*dq;	
							U(nid,isbelonged) = 0;
							U(nid,cid) = 1; 
						else
							update = true;
							U(nid,isbelonged) = 0; 
							dQ = dQ + 2*dq;	
						end
						if any(U(:,isbelonged))==false
							U(:,isbelonged) = [];
						end 
					end
				end
			end

			
			function G = second_step(G,Ct)
				G = graph(Ct'*G.adjacency_matrix()*Ct);	
			end
			
			function score = calc_node_score(C,G,q);
				A = G.adjacency_matrix();
				d = G.degree();
				n = G.sz();	

				score = sparse(n,1,0);
				for cid = 1:size(C,2)
					c = C(:,cid);
					dc = c'*d;
					score = score + (A*c/(2*m)-d.*dc/(4*m^2) +2*q).*c;
				end
				score(isnan(score))=0;
			end
	
			function P = find_periphery(G,C)
				T = G.transition_matrix();
				m = G.edgenum();
				d = G.degree();
				%P = (C'*T*diag(d/(2*m)))'-d*sum(C,1)/(2*m);% - d*sum(C,1)/(2*m);
				
				P = (C'*T)' - ones(n,1)*(d'*C)/(2*m);
				%P = (C'*T*diag(d/(2*m)))' - d*(d'*C)/(4*m^2);
				P(P<0)=0;	
				if size(C,2)==1
					P =sign(P); 
					P(any(C,2)) =0;
					return
				end
					
				[val,cid] = max(P');
				
				P = sparse(1:n,cid,1,size(T,1),size(C,2));
				P(any(C,2) | val'==0,:) =0;
				
			end
		end
		
		function model = init_model(self,model)
			if ~isfield(model,'lambda');model.lambda = 1;end
			if ~isfield(model,'itnum');model.itnum = 10;end
			if ~isfield(model,'solver');model.solver = 'louvain';end
		end
	end
end
