classdef birkan 
	methods (Access = public)
		function [C,P,Q,Qs,score,model] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
			C= [];  P = [];Q = [];score = [];
			[C,P,Q,score] = self.EMalgorithm(A,model);
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
	end
	methods (Access = private)
		function C = initializeCores(self,N)
			C = rand(N,1);
			C(C>0.5) = 1;
			C(C<=0.5) = 0;
		end
		function Z = initializeCommunities(self,N,K)
			if K==1
				Z = ones(N,1);
			else
				lab = randsample(K,N,true);
				Z = sparse(1:length(lab),lab,1,N,K);
			end
		end
		
    		function Blist = generativeCores(self, c, Z)
			Blist = cell(size(Z,2),1);
			C = c*c' + c*(1-c)' + (c*(1-c)')';
			for k = 1:size(Z,2)
				Blist{k} = (Z(:,k)*Z(:,k)').*C;
			end	
		end
		
		function [norm, logL, jlogL]=likelihood(self, A,Blist, Z,model)
			
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
        		a =full(sum(A,2)); 

        		afl = gammaln(a+1);
        		Afl = sum(gammaln(full(A+1)),2);

        		bgl = sparse(N, K);
        		abgl = sparse(N, K);
        		ABgl = sparse(N, K);
        		for k = 1:K
            			Bk = Blist{k}; 
            			b = full(sum(Bk,2));
					
            			tmp = gammaln(b);
            			tmp(b<1e-1) = 0;
			        bgl(:,k) = tmp;
			
			        AB = A+Bk;
			        ab = full(sum(AB,2));
			        tmp = gammaln(ab);
			        tmp(isnan(tmp)) = 0;
			        tmp(ab<1e-1) = 0;
			        abgl(:,k) = tmp;
			
			        tmp = gammaln(full(AB)) - gammaln(full(Bk));
			        tmp(isnan(tmp)) = 0;
			        tmp(tmp<1e-12) = 0;
			        ABgl(:,k) = sum(tmp,2);
			end
        		
			%log-likelihood of individual nodes for each community
        		logL = repmat(afl,1,K) - repmat(Afl,1,K) + bgl - abgl + ABgl + ones(N,1)*pr; 

        		%Joint log-likelihood of observations and communities
        		%see section "Inferring Meso-Scale Structures" for details
        		if ~isempty(Z)
            			norm = log(sum(logL,2)); %summation over communities. usefull in hybrid models
            			jlogL = sum(sum(Z .* logL));
        		else
            			norm = logL;
            			jlogL = sum(logL) %summation of individual log-likelihoods
			end

        		%return norm, logL, jlogL
   		end
		function c1 = updaterCores(self, c0)
		        ind = randsample(length(c0),1);
		        c1 = c0;
		        c1(ind) = 1 - c1(ind);
		
		        %When using a non-binary model (coreness ranges between 0-1) use below two lines
		        %update = np.random.normal(0, 0.1)
		        %c1[ind] = min(1., max(0., c1[ind] + update))
		
		        %return c1
		end	
		 
		function [norm, logL, f0, c0] = maximizeCores(self, A,c0, Z, model)
        		%get the initial expected number of interactions
        		Blist = self.generativeCores(c0, Z);
        		%Joint likelihood of initial observations and assignment
        		%If this is a pure core-periphery structure without any communities, Z is an array of ones
			
        		[norm, logL, f0] = self.likelihood(A,Blist, Z,model);
			counter = 0;	
        		for it = 1:model.maxIt
            			c1 = self.updaterCores(c0);
            			Blist = self.generativeCores(c1, Z);
            			[tmpnorm, tmplogL, f1] = self.likelihood(A,Blist,Z,model);
				
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
			                disp(sprintf('It: %0.4d S: %0.4d logL: %0.2f',it, counter, full(f0)))
            			elseif counter > 750
                			disp(sprintf('It: %0.4d S: %0.4d logL: %0.2f',it, counter, full(f0)))
                			break
				end
			end

        		%return norm, logL, f0, c0
		end
		
		function [C,P,Q,score] = EMalgorithm(self,A,model)
			
			Cbest = [];Zbest = [];
			bestLogL = -1e+10;
			N = size(A,1);
			for rep = 1:model.numRep
			
				c = self.initializeCores(N);
				Z = self.initializeCommunities(N,model.K);
				
            			%Solve mixture model using EM algorithm
            			f0 = 0;
            			if model.K==1; maxInIt = 1;end
				for it=1:model.maxInIt
                			%For a fixed Z optimize model parameters i.e. coreness vector
                			[norm, logL, f1, c] = self.maximizeCores(A,c, Z, model); % norm: Nx1, logL: NxK, f1:1x1
					
                			%Then update Z
                			if model.K~=1	
                    				Z = real(exp(logL - repmat(norm,1,model.K)));
							
						%denom = 1./sum(Z,2);
						%denom(isnan(denom)) = 0;denom(isinf(denom)) = 0;
						%Z = diag(denom)*Z;
					end
                			diff = abs(f1-f0);
                			f0 = f1;

                			if diff <= model.tol
                    				break
					end	
				end	
            			if f0 > bestLogL
                			bestC = c;
                			bestZ = Z;
                			bestLogL = f0;
				end
            			disp(sprintf('Repeat:%d Best log-likelihood: %0.2f', rep, bestLogL))
			end
			
			% discritize Z
			C = diag(Cbest)*Z;
			P = diag(1-C)*Z;
			Q = bestLogL;
			score = sum(Z,2);	
		end	
		function model = init_model(self,model)
			if ~isfield(model,'K');model.K=1;end
			if ~isfield(model,'numRep');model.numRep=10;end
			if ~isfield(model,'maxInIt');model.maxInIt=5000;end
			if ~isfield(model,'maxIt');model.maxIt=5000;end
			if ~isfield(model,'tol');model.tol=1e-6;end
			if ~isfield(model,'coreType');model.coreType=1;end
		end
		
	end
end
