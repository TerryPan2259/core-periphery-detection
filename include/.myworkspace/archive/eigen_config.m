%addpath /panfs/panasas01/emat/sk16722/program/lib/community
%addpath /panfs/panasas01/emat/sk16722/program/lib/graph
%addpath(genpath('/panfs/panasas01/emat/sk16722/program/lib/YALMIP'));
	
classdef eigen_config 
	methods (Access = public)
		function [C,P,Q,Qs,score,model] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model,A);
			[C,P,Q,Qs,score] = self.hard_assignment(A,model);		
		end
			
	
		function q = eval(self,A,x)
			d = sum(A,2);	
			n = size(A,1);
			m = sum(d)/2;
			An = d*d'/(2*m);
			q = trace(x'*(A-An)*x);
		end
	end
	methods (Access = private)
		
		function [C,P,Q,Qs,score,model]=hard_assignment(self,A,model)
			d = sum(A,2);	
			n = size(A,1);
			m = sum(d)/2;
			
			function [q,qs] = eval(C,S)
				qs = zeros(size(C,2),1);
				for k = 1:size(C,2)
					qs(k) = ((S.*C(:,k))'*A*(S.*C(:,k)) - sum(S(C(:,k)>0).*d(C(:,k)>0)).^2/(2*m))/(2*m);
				end
				q = sum(qs);	
			end
			
			%S = ones(n,1);S = S/sqrt(S'*S);
			S = rand(n,1);S = S/sqrt(S'*S);
			C = model.C;
			com = community();commodel.name = 'modularity';commodel.itnum = 20;
			updated = true;qbest = -Inf;
			while updated
%				commodel.C = C;
%				
%				commodel.node_weight = S;
%				c = com.detect(A,commodel);
%		
%				q = eval(c,S);
%					
%				if qbest < q
%					qbest = q;
%					C = c;
%					updated = true;
%				else
%					break
%				end
				
				x = sdpvar(n,1);
				Constraints = [x'*x ==1, x>=0];
				J = 0;
				
				J = J + (x'*A*x - sum(x.*d).^2/(2*m))/(2*m);
				%for k = 1:size(C,2)
				%	J = J + ((x.*C(:,k))'*A*(x.*C(:,k)) - sum(x(C(:,k)>0).*d(C(:,k)>0)).^2/(2*m))/(2*m);
				%end
					
				options = sdpsettings('verbose',0,'usex0',1);
				assign(x,S);
				sol = optimize(Constraints,-J,options);
				s = value(x);
					
				
				q = eval(C,s);
				if qbest < q
					qbest = q;
					S = s;
						
					updated = true;
				else
					break
				end
			end
			[Q,Qs] = eval(C,S);	
			score = S;
			C = diag(S)*C;
			C = C * sparse(1:size(C,2),1:size(C,2),1./sqrt(sum(C.^2,1)));
			P = 1-C;
		end
		
		function [C,P,Q,Qs,score,model]=als(self,A,model)
			d = sum(A,2);	
			n = size(A,1);
			m = sum(d)/2;
			An = d*d'/(2*m);
			Mod = (A-An)/(2*m);
		
			function [c,ceq] = mycon(x)
				c = 0;
				ceq = sqrt(x'*x)-1;
			end
		
		
			function qval = eval(Q,lam,X)
				tmp = Q;
				for k = 1:size(X,2)
					tmp = tmp - lam(k) * X(:,k)*X(:,k)';
				end
				qval = sum(tmp(:).^2);
			end
			
		
			C = [];Qs = [];
			K = model.K;
			X = model.X;
			lams = rand(K,1);
			%Q = A;	
			oldq = Inf;
			
			function [y,grady] = quadobj(x,Q)
				y = -1/2*x'*Q*x;
				if nargout > 1
				    grady = -Q*x;
				end
			end
			
			function [y,yeq,grady,gradyeq] = quadconstr(x,d)
				y = x'*x-1;
				yeq = [];
				%yeq = x'*x-1;
				grady = 2*x;
				%gradyeq = 2*x;
				gradyeq = [];
			end	
			
			function hess = quadhess(x,lambda,Q)
				hess = Q;	
			end	
			
			
				
			while true
				for k =1:K
					R = zeros(n,n);
					for l =1:K
						if k==l;continue;end
						R = R + lams(l)*X(:,l)*X(:,l)';	
					end
					Q = (Mod-R);
					options = optimoptions(@fmincon,'Algorithm','interior-point',...
    						'GradObj','on','GradConstr','on','Hessian','user-supplied',...
    						'HessFcn',@(x,lambda)quadhess(x,lambda,Q));
				
					[x,fval] = fmincon(@(x) (-x'*Q*x),X(:,k),[],[],[],[],zeros(n,1),ones(n,1), @mycon);
					%fun = @(x) quadobj(x,Q);
					%nonlconstr = @(x) quadconstr(x,-1);
					
					%[x,fval,eflag,output,lambda] = fmincon(fun,X(:,k),...
    					%				[],[],[],[],zeros(n,1),ones(n,1),nonlconstr,options);
					a = sqrt(x'*x);
					X(:,k)= x/a;
					lams(k) = X(:,k)'*Q*X(:,k);
				end
				q = eval(Q,lams,X);	
				if oldq <=q
					break;	
				end
				oldq = q;	
			end	
			C = X;	
			Qs = lams;
			Q = sum(Qs);
			score = sum(C,2);
			P = C;
			%[C,P,Q,score] = self.gradient(graph(A));
			%Qs = Q;
		end
		
		function model = init_model(self,model,A)
			n = size(A,1);
			if ~isfield(model,'disp');model.disp = false;end
			if ~isfield(model,'K');model.K= 1;end
			if ~isfield(model,'w');model.w= 1/sum(A(:));end
			if ~isfield(model,'type');model.type='hard';end
			if ~isfield(model,'parallel');model.parallel=false;end
			
			if ~isfield(model,'C');
				cid = randsample(1:model.K,size(A,1),true);
				model.C  = sparse(1:length(cid),cid,1);
			end
	sparse(1:n,1:n,1);
			if ~isfield(model,'X');
				model.X= rand(n,model.K);
				model.X= model.X*diag(1./sqrt(sum(model.X.^2)));
			else
				if size(model.X,2) <model.K
					tmp = rand(n,model.K-size(model.X,2));
					tmp= tmp*diag(1./sqrt(sum(tmp.^2)));
					model.X = [model.X,tmp]; 
				end
			end
		end
	end
end
