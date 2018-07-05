classdef birkan_cp 
	methods (Access = public)
		function [C,P,Q,Qs,score,model,dt] = detect(varargin)
			self = varargin{1};
			G = varargin{2};
			A = G.adjacency_matrix();
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.initParam(model);
			C= [];  P = [];Q = [];score = [];
			Q = -Inf;
			dt = cputime;
			for rep = 1:model.numRun
			for k = model.Kmin:model.Kmax
				model.K = k;
				[Ct,Pt,Qt,scoret] = self.EMalgorithm(A,model);
				if Qt > Q
					C = Ct;P = Pt;
					Q = Qt;	
					score = scoret;
				end 
			end
			end	
			slice = sum(C + P,1)>0;
			C = C(:,slice); P = P(:,slice);
			
			dt = cputime-dt;
			Qs = Q;
		end
	
		function q = eval(self,A,x)
			n = size(A,1);
			L = nchoosek(1:n,2);
			eids = sub2ind([n,n],L(:,1),L(:,2));
			a = A(eids);	
			xx = x(L(:,1)) + x(L(:,2)) - ( x(L(:,1)) .* x(L(:,2)) );	
			q = corr(a,xx);	
			
			if isnan(q)
				q = 0;
			end	
		end
	
		function model = initParam(self,model)
			if ~isfield(model,'K');model.K=2;end
			if ~isfield(model,'Kmin');model.Kmin=1;end
			if ~isfield(model,'Kmax');model.Kmax=5;end
			if ~isfield(model,'a');model.a=1;end
			if ~isfield(model,'b');model.b=0;end
			if ~isfield(model,'maxInIt');model.maxInIt=50;end
			if ~isfield(model,'maxIt');model.maxIt=100;end
			if ~isfield(model,'tol');model.tol=1e-6;end
			if ~isfield(model,'coreType');model.coreType=1;end
			if ~isfield(model,'disp');model.disp=false;end
		end
	end
	methods (Access = private)
		function [C,P,Q,score,dt] = original(self,A,model)
			
			str = python('birkan.py',mat2str(full(A)),num2str(full(model.K)),num2str(full(model.itnum)));
			R = JSON.parse(str);
				
			C = cell2mat(R.C)';
			Z = cell2mat(cellfun(@(x) [x{1,:}], R.Z(:),'UniformOutput',false));
			Q = R.Q;
			dt = R.dt;
			% discritize Z
			if size(Z,2)==1
				U = ones(size(Z,1),1);	
			else
				U = sparse(size(Z,1),size(Z,2),0);
				for i = 1:size(Z,1)
					[~,idx] = max(Z(i,:));
					U(i,idx) = 1;
				end
			end
			
			C = diag(C)*U;
			P = U-C;
			score = sum(Z,2);	
		end
	
		function C = initializeCores(self,A,N)
			deg = sum(A,2);
			C = rand(N,1);
			C(deg>mean(deg)) = 1;
			C(deg<=mean(deg)) = 0;
		end
		function Z = initializeCommunities(self,N,K)
			if K==1
				Z = ones(N,1);
			else
				lab = randsample(K,N,true);
				Z = sparse(1:length(lab),lab,1,N,K);
			end
		end
		
    		function BB = generativeCores(self, cc, ZZ,a,b)
			BB = cell(size(ZZ,2),1);
			C = cc*cc' + cc*(1-cc)' + (cc*(1-cc)')';
			for k = 1:size(ZZ,2)
				BB{k} = a*(ZZ(:,k)*ZZ(:,k)').*C+b;
			end
		end
		function [norm,logL,jlogL] = likelihood(self,A,B,Z,model)
			K = model.K;
			N = size(A,1);
        		if ~isempty(Z) %if there are communities
            			K = model.K;
            			pr = sum(Z,1);
            			pr = pr / sum(pr); %prior probabilities of communities
        		else %if not
            			K = 1;
            			pr = [1];
			end
        		%This part calculates the Dirichlet-Multinomial log-likelihood (see equation 3 of the paper)
        		afl = abs(gammaln(sum(A,2)+1));
				
        		%Afl = sum(abs(gammaln(full(A+1))),2); % I comment this out since Afl is always a zero matrix for unweighted networks 
        		bgl = sparse(N, K);
        		abgl = sparse(N, K);
        		ABgl = sparse(N, K);
        		for k = 1:K
            			Bk = B{k}; 
            			b = full(sum(Bk,2));
					
            			tmp = abs(gammaln(b));
            			tmp(b<1e-1) = 0;
			        bgl(:,k) = tmp;
			
			        AB = A+Bk;
			        ab = full(sum(AB,2));
			        tmp = abs(gammaln(ab));
			        tmp(isnan(tmp)) = 0;
			        tmp(ab<1e-1) = 0;
			        abgl(:,k) = tmp;
			
			        tmp = abs(gammaln(full(AB))) - abs(gammaln(full(Bk)));
			        tmp(isnan(tmp)) = 0;
			        tmp(tmp<1e-12) = 0;
			        ABgl(:,k) = sum(tmp,2);
			end
        		
			%log-likelihood of individual nodes for each community
        		logL = repmat(afl,1,K) + bgl - abgl + ABgl + log(ones(N,1)*pr); 
        		%logL = repmat(afl,1,K) - repmat(Afl,1,K) + bgl - abgl + ABgl + log(ones(N,1)*pr); 
			
        		%Joint log-likelihood of observations and communities
        		%see section "Inferring Meso-Scale Structures" for details
        		if ~isempty(Z)
            			norm = self.logsum(logL); %summation over communities. usefull in hybrid models

            			jlogL = sum(sum(Z .* logL));
        		else
            			norm = logL;
            			jlogL = sum(logL) %summation of individual log-likelihoods
			end
		end	
%		function [norm, logL, jlogL]=likelihood(self, A,deg,Blist, Z,model)
%		%function [norm, logL, jlogL]=likelihood(self, A,deg,Blist, Z,model)
%			
%			K = model.K;
%			N = size(A,1);
%        		if ~isempty(Z) %if there are communities
%            			K = model.K;
%            			pr = sum(Z,1);
%            			pr = pr / sum(pr); %prior probabilities of communities
%        		else %if not
%            			K = 1;
%            			pr = [1];
%			end
%        		%This part calculates the Dirichlet-Multinomial log-likelihood (see equation 3 of the paper)
%        		afl = abs(gammaln(deg+1));
%				
%        		%Afl = sum(abs(gammaln(full(A+1))),2); % remove comment out if A is weighted matrix
%			Afl = sum(gammaln(full(A+1)),2);
%        		bgl = sparse(N, K);
%        		abgl = sparse(N, K);
%        		ABgl = sparse(N, K);
%        		for k = 1:K
%            			Bk = Blist{k}; 
%            			b = full(sum(Bk,2));
%					
%            			tmp = abs(gammaln(b));
%            			tmp(b<1e-1) = 0;
%			        bgl(:,k) = tmp;
%			
%			        AB = A+Bk;
%			        ab = full(sum(AB,2));
%			        tmp = abs(gammaln(ab));
%			        tmp(isnan(tmp)) = 0;
%			        tmp(ab<1e-1) = 0;
%			        abgl(:,k) = tmp;
%			
%			        tmp = abs(gammaln(full(AB))) - abs(gammaln(full(Bk)));
%			        tmp(isnan(tmp)) = 0;
%			        tmp(tmp<1e-12) = 0;
%			        ABgl(:,k) = sum(tmp,2);
%			end
%        		
%			%log-likelihood of individual nodes for each community
%        		logL = repmat(afl,1,K) - repmat(Afl,1,K) + bgl - abgl + ABgl + log(ones(N,1)*pr); 
%			
%        		%Joint log-likelihood of observations and communities
%        		%see section "Inferring Meso-Scale Structures" for details
%        		if ~isempty(Z)
%            			norm = self.logsum(logL); %summation over communities. usefull in hybrid models
%
%            			jlogL = sum(sum(Z .* logL));
%        		else
%            			norm = logL;
%            			jlogL = sum(logL) %summation of individual log-likelihoods
%			end
%
%   		end
		
		function s = logsum(self,array)
		   amax = max(array')';
		   
		   s = log(sum(exp(array - amax*ones(1,size(array,2))),2)) + amax;
		
		end
		
		function cc1 = updaterCores(self, cc0)
		        ind = randsample(length(cc0),1);
		        cc1 = cc0;
		        cc1(ind) = 1 - cc1(ind);
		
		        %When using a non-binary model (coreness ranges between 0-1) use below two lines
		        %update = np.random.normal(0, 0.1)
		        %c1[ind] = min(1., max(0., c1[ind] + update))
		
		        %return c1
		end	
		 
		function [norm, logL, f0, c0] = maximizeCores(self, A,c0, Z, model)
        	
			%get the initial expected number of interactions
        		B = self.generativeCores(c0, Z,model.a,model.b);
        		[norm, logL, f0] = self.likelihood(A,B,Z,model);
			%Blist = self.generativeCores(c0, Z,model.a,model.b);
		
		
        		%Joint likelihood of initial observations and assignment
        		%If this is a pure core-periphery structure without any communities, Z is an array of ones
			counter = 0;	
        		for it = 1:model.maxIt
            			c1 = self.updaterCores(c0);
            			B = self.generativeCores(c1, Z,model.a,model.b);
            			%Blist = self.generativeCores(c1, Z,model.a,model.b);
            			[tmpnorm, tmplogL, f1] = self.likelihood(A,B, Z,model);
				
            			%This is a very primitive deterministic optimization
            			%Consider using a better method such as Simulated Annealing
            			%by using a random process to calculate 'test'
            			test = f1 > f0;
            			counter = counter + 1;
            			if test
                			c0 = c1;
			                f0 = f1;
					norm = tmpnorm;
					logL = tmplogL;
			                counter = 0;
			                if model.disp;disp(sprintf('It: %0.4d S: %0.4d logL: %0.2f',it, counter, full(f0)));end
            			elseif counter > 750
                			if model.disp;disp(sprintf('It: %0.4d S: %0.4d logL: %0.2f',it, counter, full(f0)));end
                			break
				end
			end
        		%return norm, logL, f0, c0
		end
		
		function [C,P,Q,score,dt] = EMalgorithm(self,A,model)
			ts = cputime;	
			Cbest = [];Zbest = [];
			bestLogL = -1e+10;
			N = size(A,1);
			deg = full(sum(A,2));
			c = self.initializeCores(A,N);
			Z = self.initializeCommunities(N,model.K);
			
            		%Solve mixture model using EM algorithm
            		f0 = 0;
			maxInIt = model.maxInIt;
            		if model.K==1; maxInIt = 1;end
			for it=1:maxInIt	
               			%For a fixed Z optimize model parameters i.e. coreness vector
               			[norm, logL, f1, c] = self.maximizeCores(A,c, Z, model); % norm: Nx1, logL: NxK, f1:1x1
               			%Then update Z
               			if model.K~=1	
                   				Z = exp(logL - repmat(norm,1,model.K));
				end
					
               			df = abs(f1-f0); % I do not know why we have to take abs() for the value; It should be always positive.
               			f0 = f1;
               			if df <= model.tol
                   				break
				end	
			end	
		
            		if f0 > bestLogL
                		bestC = c;
                		bestZ = Z;
                		bestLogL = f0;
			end
			
			% discritize Z
			if size(Z,2)==1
				U = ones(size(Z,1),1);	
			else
				U = sparse(size(bestZ,1),size(bestZ,2),0);
				for i = 1:size(bestZ,1)
					[~,idx] = max(bestZ(i,:));
					U(i,idx) = 1;
				end
			end
			C = diag(bestC)*U;
			P = diag(1-bestC)*U;
			Q = bestLogL;
			score = sum(Z,2);	
			dt = cputime-ts;	
		end	
		
	end
end
