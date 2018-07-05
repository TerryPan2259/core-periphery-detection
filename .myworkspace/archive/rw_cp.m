classdef rw_cp
	methods (Access = public)
		function [C,P,Q,Qs,score,model] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
		
			
			C = []; P = [];Qs = [];score = []; Q = -Inf;
			for it = 1:model.itnum	
				switch model.solver
					case 'label_switching'
						[Ct,Pt,Qt,Qst,scoret] = self.label_switching(graph(A),model);
					otherwise	
						disp('unknown solver for rwcp')
				end
				
				if Qt > Q
					C = Ct;P = Pt;Q = Qt;Qs = Qst;score = scoret; 
				end
			end
		end
	end
	methods (Access = private)
		function [C,P,Q,Qs,score] = label_switching(self,G,model)
			
			% ******************************************* 
			% Main Routine 
			% ******************************************* 
			
			% --------------------------------
			% Initialise  
			% --------------------------------
			T = G.transition_matrix();
			n = G.sz();
			
			% --------------------------------
			% Label switching (Main routine) 
			% --------------------------------
			C = sparse(1:n,1:n,1);P = sparse(n,n,0);
			dQ = 0; 
			ord = randsample(n,n)';
			updated = true;
			
			while(updated)
				updated = false;
				for nid =ord
					iscore = any(C(nid,:));
					Nc = sum(C,1);
					Np = sum(P,1);
					cQv = (C'*T(:,nid))' - model.lambda * (Nc-C(nid,:))/(n-1);	
					pQv = (P'*T(:,nid))' - model.lambda * (Np-P(nid,:))/(n);	
					vQc = T(nid,:)*C - model.lambda * (Nc-C(nid,:))/(n-1);
					vQp = T(nid,:)*P - model.lambda * (Np-P(nid,:))/(n-1);
					vQv = 0;
					%vQv = T(nid,nid)-1/n;
					
					% core join gain 
					%cjgain = T(nid,:)*(C+P) - 1 + T(nid,nid);
					cjgain = cQv + pQv + vQc + vQp + vQv;
					
					% periphery join gain 
					pjgain = cQv + pQv + vQc + vQp -vQp - pQv;
					%pjgain = (C'*T(:,nid))';
					% calculate diff of Q vals
					[v2c,v2cid] = max(cjgain);
					[v2p,v2pid] = max(pjgain);
					
					
					if iscore
						% leave gain
						cid = find(C(nid,:)>0);
						lgain = -cQv - pQv - vQc - vQp + vQv;
						lgain = lgain(cid);
					else
						% leave gain
						cid = find(P(nid,:)>0);
						lgain = -cQv - pQv - vQc - vQp + vQp + pQv;
						lgain = lgain(cid);
					end
					
					if v2p < v2c & eps*1000 < v2c + lgain & cid ~=v2cid
						C(nid,:) = 0; 
						P(nid,:) = 0; 
						C(nid,v2cid) = 1; 
						updated = true;
					elseif v2c <= v2p & eps*1000 < v2p + lgain & cid ~=v2pid
						C(nid,:) = 0; 
						P(nid,:) = 0; 
						P(nid,v2pid) = 1; 
						updated = true;
					end
				end	
			end
			slice = any(C+P,1);
			C = C(:,slice);
			P = P(:,slice);
			
			Nc = sum(C,1); Np = sum(P,1);
			score = (sum( (T*(C+P)).*C,2) - C*(Nc+Np)'/n)/n;
			Qs = (diag((C+P)'*T*(C+P))' - ((Np+Nc).^2 - Np.^2 -Nc)/n)/n;
			Q = sum(Qs);
			 
		end
		
		function model = init_model(self,model)
			if ~isfield(model,'lambda');model.lambda = 1;end
			if ~isfield(model,'itnum');model.itnum = 10;end
			if ~isfield(model,'solver');model.solver = 'label_switching';end
		end
	end
end
