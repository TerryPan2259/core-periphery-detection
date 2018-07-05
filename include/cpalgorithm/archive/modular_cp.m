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
			q = model.lambda * (sum(T,2) - ones(n,1));
			
			%q = model.lambda * (sum(T,2) - ones(n,1))/n;
			%q = model.lambda * (sum(T,2) - d*n/(2*m))/n;
			
			%q = model.lambda * (sum(T,2) - ones(n,1))/n;
			q(~any(T,1))=0;
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
				break	
				%Gt = second_step(Gt,Ct);
			end
			
			if size(C,2)==0;
				C = sparse(G.sz(),1,0);
				P = sparse(G.sz(),1,0);
				score = sparse(G.sz(),1,0);
				Qs = [0]; Q = 0;
				return;
			end
			P = find_periphery(T,C,G.edgenum(),G.degree());
			U = C + P;
			
			com = cpdetection();
			mymodel.name = 'borgatti';
			mymodel.solver = 'eigen';
			C = [];P = [];
			n = size(A,1);
			Ag = G.adjacency_matrix();	
			for i = 1:size(U,2)
				u = U(:,i);
				sz = sum(u);
				As = Ag(u>0,u>0);
				R = sparse(find(u),1:sz,1,n,sz);
				
				[D,dQ]= com.detect(As,mymodel);
				if length(D)==0
					C = [C,R*ones(size(As,1),1)]; P = [P,zeros(size(A,1),1)]; 
				else
					C = [C,R*D]; P = [P,R*(1-D)]; 
				end
			end
			
			score = calc_node_score(C,G,q);
			Qs = score'*C;Q = sum(Qs);
			[Qs,ord] = sort(Qs,'descend');
			
			%P(:,sum(C,1)<=1)=[];
			%C(:,sum(C,1)<=1)=[];

	
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
					
				M = (A - G.density()*(ones(n)-eye(n)));
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
							qr = -qp(isbelonged) - 2*q(nid) - vMv(nid);
							
							%qr = -q(nid) - ( (vAc(isbelonged) + vAv(nid) - d(nid)*(d(nid)+dU(isbelonged))/(2*m))/(2*m)  ); 
							%qp = vAc/(2*m) -d(nid)*dU/(4*m^2);
			
							qp = qp - qp(isbelonged);	
						else
							vAc = A(nid,:)*U;
							dU = d'*U;
							qr = 0;
								
							qp = M(nid,:)*U + 2*q(nid) + vMv(nid);
								
							%qp = vAc/(2*m) -d(nid)*dU/(4*m^2) + q(nid);
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
				m = G.edgenum();
				score = sparse(n,1,0);
				if m==0;return;end;
				for cid = 1:size(C,2)
					c = C(:,cid);
					dc = c'*d;
					score = score + (A*c/(2*m)-d.*dc/(4*m^2) +2*q).*c;
				end
			end	
			function P = find_periphery(T,C,m,d)
				if m==0;P = sparse(size(T,1),1,0);return;end
				P = (C'*T)';% - d*sum(C,1)/(2*m);
				%P = (C'*T)'-ones(n,1)*(C'*d)'/(2*m);% - d*sum(C,1)/(2*m);
				%P = (C'*T*diag(d/(2*m)))'-d*sum(C,1)/(2*m);% - d*sum(C,1)/(2*m);
				%P = (C'*T)'/n - ones(n,1)*(sum(T,2)'*C)/(n*n);
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
