classdef cpdetection
	methods (Access = public)
		 function [C,P,Q,Qs,score,model] = detect(varargin)
			self = varargin{1};
			A = varargin{2};
			model = [];
			if nargin >2
				model = varargin{3};
			end
			model = self.init_model(model);
			
			switch model.name
				case 'borgatti'
					cpd = borgatti_cp();
				case 'borgatti_diag'
					cpd = borgatti_cp_diag();
				case 'borgatti_multi'
					cpd = borgatti_multi();
				case 'eigen'
					cpd = eigen_cp();				
				case 'eigen_config'
					cpd = eigen_config();
				case 'rombach'
					cpd = rombach_cp();
				case 'rombach_multi'
					cpd = rombach_multi();
				case 'richcore'
					cpd = richcore();				
				case 'fabio'
					cpd = fabio_cp();	
				case 'recursive'
					cpd = recursive_cp();% appy single core/periphery detection recursively	
				case 'com_cp' 
					cpd = community_cp();% divide community and then divide core and periphery	
				case 'modular'
					cpd = modular_cp();
				case 'rw_cp'
					cpd = rw_cp();
				case 'hook_potts_cp'
					cpd = hook_potts_cp();
				otherwise
					display ('unknown model name in cpdetection')
			end
			[C,P,Q,Qs,score,model] = cpd.detect(A,model);
			if isfield(model,'pval')
				slice = self.permutation_test(C,P,A,model);
				C = C(:,slice);P = P(:,slice);Qs = Qs(slice);Q = sum(Qs);
			end
			%if size(C,2)>1	
			%	remove = any(C,1)==false;
			%	C(:,remove) =[];P(:,remove) =[];Qs(remove) = [];Q= sum(Qs);
			%end	
		end
		
		function pval = calc_pval(varargin)
			C = varargin{2};
			As = varargin{3};
			sampleNum = 1000;
			if nargin ==4
				sampleNum = varargin{4};
			end
			cpmodel.name = 'kl';	
			cpd = borgatti_cp();
			
			
			% preparation	
			n = size(As,1);L = nchoosek(1:n,2);eids = sub2ind([n,n],L(:,1),L(:,2));
			a = As(eids);	
			q = cpd.eval(As,C);
			count = 0;
			for t = 1:sampleNum % pretest
				At = sparse(L(:,1),L(:,2),a(randsample(length(a),length(a))),n,n);% random permutation 
				At = At + At'; % symmetrise
				[~,~,qt]=cpd.detect(At,cpmodel);
				if q <= qt
					count = count + 1;
				end
			end
			pval = count / sampleNum;
		end
		
		function significant = permutation_test(self,C,P,A,model)
			if model.disp
				disp('start statistical test...');
			end
			if isfield(model, 'sampleNum')
				sampleNum = model.sampleNum;
			else
				sampleNum = 3000;	
			end
			clear cpmodel;		
			cpmodel.name = 'kl';	
			cpd = borgatti_cp();
			significant = true(1,size(C,2));
		
			corrected_pval = 1-(1-model.pval)^(1/size(C,2));
			
			for cid = 1:size(C,2)
				if sum(C(:,cid)+P(:,cid),1)<=2;significant(cid)=false;continue;end
					
				slice = any(C(:,cid)+P(:,cid),2);
				As = A(slice,slice); % slice
				
				% preparation	
				n = size(As,1);L = nchoosek(1:n,2);eids = sub2ind([n,n],L(:,1),L(:,2));
					
				a = As(eids);	
				
				q = cpd.eval(As,C(slice,cid));
				count = 0;
				if model.disp
					disp(sprintf('testing .... n=%d, m=%d',n,full(sum(As(:))/2)));
				end
				Qt = zeros(sampleNum,1);
				for t = 1:sampleNum % pretest
					At = sparse(L(:,1),L(:,2),a(randsample(length(a),length(a))),n,n);% random permutation 
					At = At + At'; % symmetrise
					[~,~,qt]=cpd.detect(At,cpmodel);
					if q <= qt
						count = count + 1;
						
						if count/sampleNum >corrected_pval
							significant(cid) = false;
							break	
						end
					end
				end
				
				
				if model.disp
					disp(sprintf('%d/%d stest finished (res=%d)',cid,size(C,2),significant(cid)))
				end
			end	
		end
	end
	
	methods (Access = private)
		function model = init_model(self,model)
			if ~isfield(model,'name');model.name = 'borgatti';end
			if ~isfield(model,'disp');model.disp = false;end
			if ~isfield(model,'palarell');model.palarell = false;end
		end
		
		function Alist = ergraph(self,rho,N,T)
			M = N*(N-1)/2;
			L = nchoosek(1:N,2);
			t = 1;
			Alist = cell(T,1);	
			while(t <=T)	
				m = binornd(M,full(rho));
				eidx = randsample(M,m,false);
				A = sparse(L(eidx,1),L(eidx,2),1,N,N);
				A = A + A'; 
				Alist{t,1} = A;
				t = t + 1;
			end
		end
	
	end
end	
