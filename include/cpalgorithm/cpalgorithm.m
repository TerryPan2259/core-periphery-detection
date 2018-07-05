% Container class for be_cp, two_step_cp and km_cp classes
% 
% @version 0.1.0
% @date 2017/02/19
%
% @author Sadamori Kojaku 
classdef cpalgorithm

	properties
		algorithms;
	end

	methods ( Access = public )

		function obj = cpalgorithm()	
			obj.algorithms = { % List of available algorithms
				% Borgatti-Everett algorithm
				'random', random_cp();
				
				% Borgatti-Everett algorithm
				'be', be_cp();
				
				% Borgatti-Everett algorithm (just a clique. pattern 2)
				'be_without_cp_edges', be_cp();
				
				% Borgatti-Everett algorithm (core-periphery edges are ignored)
				'be_continuous', be_cp();
				
				% Borgatti-Everett algorithm (core-periphery edges are ignored)
				'be_cont_dcor', be_cp();
				
 				% Two-step algorithm (detect cp and com individually, and then merge)
				'two-step', two_step_cp(); 
 				
 				% Two-step algorithm (detect coms first, then divide each of them into a core and a periphery part)
				'divisive', two_step_cp(); 
				
				% Kojaku-Masuda algorithm
				'km', km_cp();
				
				% test
				'km_cont', km_cont_cp();
				
				% Rombach algorithm 
				'rombach', rombach_cp();
				
				% Rombach algorithm (degree corrected) 
				'rombach_dcor', rombach_cp();
				
				% Fabio's algorithm 
				'dellarossa', dellarossa_cp();
				
				% discrete variant of the MINRES algorithm 
				'minres', minres_cp();
				
				% discrete variant of the MINRES algorithm 
				'minres_discrete', minres_cp();
				
				% continuous variant of the MINRES algorithm 
				'minres_cont', minres_cp();
				
				% Rich-cores  
				'richcore', richcore_cp();
				
				% Zhang-Martin algorithm 
				'zm', zm_cp();
				
				% Xian's algorithm 
				'xian', xian_cp();
				
				% Cucuringu  
				'low_rank_core', cucuringu_cp();
				
				% Cucuringu  
				'lap_core', cucuringu_cp();
				
				% Cucuringu  
				'lapsgn_core', cucuringu_cp();
				
				% Star 
				'star', star_cp();
				
				% Maxcut 
				'maxcut', maxcut_cp();
				
				% tunc-velma 
				'birkan', birkan_cp();
				
				% newman modularity (this only gives communities) 
				'modularity', modularity_wrapper();
				
				% newman modularity (this only gives communities) 
				'yan', yan_cp();
			};	
		end
		
		function [C, P, Q, Qs, score, param, cpu_time] = detect( self, A, param )
			% --------------------------------
			% Initialise
			% --------------------------------
			active = find( any(A,1) | any(A,2)');
			H = sparse(1:length(active),active,1,length(active),size(A,1));
			G = graph( H*A*H' );

			if length(active)==0
				C = [];P = [];
				Q = [];Qs = [];
				score = [];
				cputime = 0;
				return;
			end				

			% --------------------------------
			% Detect core-periphery pairs 
			% --------------------------------
			cpd = self.select_algorithm( param.name );	
			[C, P, Q, Qs, score, param, cpu_time] = cpd.detect( G, param );
			
%			% --------------------------------
%			% Estimate the significance of the detected core-periphery pairs 
%			% --------------------------------
%			if isfield( param, 'pval' )
%				[slice,pvals] = self.statistical_test( C, P, H*A*H', param );
%				if ~any(slice);
%					C = sparse( G.numNode(), 1, 0 ); P=C;
%				else
%					C = C( :, slice ); P = P( :, slice );
%				end
%				[Q, Qs, score] = cpd.eval( G, C, P, param );
%			end
			C = H'*C;
			P = H'*P;
			score = H'*score;
		end
	
		function [q, qs, score] = eval( self, A, C, P, param )
			if isempty(C) | isempty(P)
				q =[];qs = [];score = [];
				return;
			end
			cpd = self.select_algorithm( param.name );
			[q, qs, score] = cpd.eval( graph(A), C, P, param );
		end
		
		function [significant,pvals,randGraphs,Crand,Prand, param] = statistical_test( self, C, P, A, param )
			
			if sum(C(:)+P(:))==0
				significant = [];	
				randGraphs = [];
				Crand =[];
				Prand = [];
				pvals = [];
				return 	
			end
			
			if ~isfield(param, 'stat_test')
				significant = boyd_stat_test(C, P, graph(A), param);
				randGraphs = [];
				pvals = [];
				Crand = [];
				Prand = [];
				return;
			end
			
			if strcmp(param.stat_test,'none')
				significant = true(1,size(C,2));
				pvals = zeros(1,size(C,2));
				randGraphs = [];
				Crand = [];
				Prand = [];
				return;
			end
		
			
			if strcmp(param.stat_test,'boyd')
				significant = self.boyd_stat_test(C, P, graph(A), param);
				randGraphs = [];
				Crand = [];
				Prand = [];
				pvals = 1-significant;

			elseif strcmp(param.stat_test,'qtest') 
				if strcmp(param.name,'km') & strcmp(param.null_model,'config')
					
					pval = 0.05;	
					pval
					N = size(A, 1);
					[r,c] = find(C+P);
					cids = full( sparse(1,r,c));
					x = full( double(any(C,2)));
			
					[r,c,v] = find(triu(A));
					L = [r,c,v];
					[pvals, nhat, qhat] = km_config_calc_pvalue(L, cids, x, N, size(L,1), param.numRun, param.numRandGraph);
					size(pvals)	
					cpval = 1 - (1-pval)^(1/size(C,2));
					significant = pvals<= cpval;
					
					randGraphs = [];
					Crand = [];
					Prand = [];
					param.nhat = nhat;	
					param.qhat = qhat;	
				else 
					[significant,pvals,randGraphs,Crand,Prand] =self.q_test(C, P, graph(A), param);
				end
			end
				
		end
		
		% Remove this after recomputing p-values !!	
		function [significant,pvals,nref,qref] = recalc_pvalues( self, C, P, A, param, Crand, Prand, randGraphs )
			if param.disp
				disp( 'starting statistical test...' );
			end
			if isfield( param, 'pval' )
				pval = param.pval;
			else
				pval = 0.05;	
			end
			
			significant = true( 1, size( C, 2 ) );
			corrected_pval = 1 - ( 1 - pval ) ^( 1 / size( C, 2 ) ); 
			Szr = [];
			Qsr = [];
			for rid  = 1:length(randGraphs)
				Ar = randGraphs{rid};	
				Cr = Crand{rid};Pr = Prand{rid};
				szr = sum(Cr+Pr,1);
				[~,qs] = self.eval( Ar, Cr, Pr, param );
				Szr = [Szr,szr];
				Qsr = [Qsr,qs];
			end
			[~,Qs] = self.eval(A, C, P, param);
			Sz = sum(C+P,1);
			pvals = self.calc_pvals_kde(Szr',Qsr',Sz',Qs');
			significant = pvals<=corrected_pval;
			nref = Szr;
			qref = Qsr;
		end
		
		function param = init( self, name )
			cpa = self.select_algorithm( name );
			clear alg; 
			param.name =name; 
			param = cpa.initParam( param );
				
			if ~isfield( param, 'disp' ); param.disp = false; end
			if ~isfield( param, 'stat_test' ); param.stat_test = 'qtest'; end
		end
	end
	
	methods ( Access = private )

		function cpd = select_algorithm( self, name )
			cpd = [];	
			for i = 1:size( self.algorithms, 1 )
				if strcmp( self.algorithms{i, 1}, name )
					cpd = self.algorithms{i, 2};	
					break
				end
			end
		end	
	
		function nids = shuffle_order( self, ids )
			[val, ord] = sort( rand( length( ids ), 1 ) );	
			nids = ids( ord );
		end	
		
		% q-value test (koujaku)	
		function [significant,pvals,randGraphs,Crand,Prand] = q_test( self, C, P, G, param )
			% --------------------------------
			% Initialise
			% --------------------------------
			if isfield( param, 'numRandGraph' )
				numRandGraph = param.numRandGraph;
			else
				numRandGraph = 500;	
			end
			if isfield( param, 'pval' )
				pval = param.pval;
			else
				pval = 0.05;	
			end
			
			rg = randgraph();
			if isfield( param, 'rgparam' )
				rgparam = param.rgparam;	
			else
				rgparam = rg.init('config');
				rgparam.type = 'soft';
			end
			rgparam.numGraph = numRandGraph;
			
			significant = true( 1, size( C, 2 ) );
			corrected_pval = 1 - ( 1 - pval ) ^( 1 / size( C, 2 ) ); 
			A = G.adjacency_matrix();
			
			randGraphs = rg.generate( A, rgparam );
			cpparam2 = param;
			if isfield( cpparam2, 'pval' )
				cpparam2 = rmfield(cpparam2,'pval');
			end
			cpparam2.disp = false;
		
			Szr = [];
			Qsr = [];
			Crand = cell(0);
			Prand = cell(0);
			for rid  = 1:length(randGraphs)
				if param.disp
					disp( sprintf('sampling core-periphery pairs %d/%d',rid,length(randGraphs)) );
				end
				Ar = randGraphs{rid};
				[Cr,Pr,~,qs] = self.detect( Ar, cpparam2);
				Crand{rid} = Cr;Prand{rid} = Pr;
				szr = sum(Cr+Pr,1);
				Szr = [Szr,szr];
				Qsr = [Qsr,qs];
			end
			Sz = sum(C+P,1);
			[~,Qs] = self.eval( A, C, P, param);
			pvals = self.calc_pvals_kde(full(Szr'),full(Qsr'),full(Sz'),full(Qs'));
			if isempty(pvals)
				significant = [];
			else
				significant = pvals<=corrected_pval;
			end
		end
	
		
		function pvals = calc_pvals_kde(self,nref,qref,ntest,qtest)
	
			if isempty(ntest) | isempty(qtest) | isempty(nref) | isempty(qref) 
				pvals = [];
				return;
			end
	
			Sig = real( cov([nref,qref]) ) + eye(2)*1e-10;
		%Sig(1,2) = 0;Sig(2,1) = 0; % test
			mean(nref)
			mean(qref)
			sqrt(Sig(1,1))	
			sqrt(Sig(2,2))	
			sqrt(Sig(1,2))	
			h = max( length(qref)^(-1/6), eps);
			pvals = ones(1,length(qtest));
			for i = 1:length(qtest)
				
				qbar = qref + Sig(1,2)/Sig(1,1) * (  ntest(i) - nref );
				
				cum = normcdf(  full( sqrt(Sig(1,1)) * (qtest(i) - qbar )/( sqrt(Sig(1,1)*Sig(2,2)-Sig(1,2)^2)*h)) );
				w = - (ntest(i)-nref).^2 /(2*h^2*Sig(1,1));
				
				w = exp( w ); % shift w to prevent overflows exp(w)
				%w = exp( w ); % shift w to prevent overflows exp(w)
					
				w = w + eps;% in case of any(w)==false 
				pvals(i) = 1-sum(cum.*w)/sum(w);
			end
		end	
		
		% Boyd 2006 Soc. Netw.		
		function [significant,pvals] = boyd_stat_test( self, C, P, G, param )
			
			% --------------------------------
			% Initialise
			% --------------------------------
			if param.disp
				disp( 'starting statistical test...' );
			end
			if isfield( param, 'numRandGraph' )
				numRandGraph = param.numRandGraph;
			else
				numRandGraph = 1000;	
			end
			if isfield( param, 'pval' )
				pval = param.pval;
			else
				pval = 0.05;	
			end
			significant = true( 1, size( C, 2 ) );
			pvals = ones( 1, size( C, 2 ) );
			
			corrected_pval = 1 - ( 1 - pval ) ^( 1 / size( C, 2 ) ); 
			 
			% --------------------------------
			% Estimate the significance of each core-periphery pair 
			% --------------------------------
			A = G.adjacency_matrix('binary');
			A = A - diag(diag(A));
			
			be = self.select_algorithm( 'be' );
			beparam = self.init('be');		
			for k = 1:size( C, 2 ) 
				if sum( C( :, k ) + P( :, k ), 1 )<=2; significant( k )=false; continue; end
				nids = any( C( :, k ) + P( :, k ), 2 ); 
				As = A( nids, nids ); 
				Qs = be.eval( graph(As), C( nids, k ), 1-C(nids,k), beparam );
				
				if param.disp
					disp( sprintf( 'testing core-periphery pair %d', k ) );
				end
				
				[r,c,v] = find(triu(As,1));	
			
				%pval = boyd_stest( [r,c], size(As,1), length(r), C(nids,k), corrected_pval, numRandGraph );
				pval = self.estpval_boyd_stest( As, Qs, corrected_pval, numRandGraph );
				significant( k ) = pval < corrected_pval;
				pvals(k) = pval;
			end
		end
	
		function estpval = estpval_boyd_stest( self, A, Q, pval, numRandGraph )
			rg = randgraph();
			rgparam = rg.init('erdos-renyi');
			rgparam.type = 'hard';
			rgparam.numGraph = 1;
			
			p = nnz(A)/ (size(A,1) * (size(A,1)-1));
			if p > 0.5 % for dense 
				A = 1-A;
				A = A - diag(diag(A));
				
				testparam = self.init('be_without_cp_edges');
				testparam.numRun = 1;
				test_alg = self.select_algorithm( testparam.name );
			else % for sparse
				testparam = self.init('be');
				testparam.numRun = 1;
				test_alg = self.select_algorithm( testparam.name );
			end
			counter = 0;
			for t = 1:numRandGraph	
				At = rg.generate( A, rgparam );
				[Ct,~,Qt] = test_alg.detect( graph(At), testparam );
				if Q <= Qt	
					counter = counter + 1;	
					if counter/numRandGraph > pval
						estpval = counter / t;
						return;
					end 
				end
			end
			estpval = counter / numRandGraph;
		end	
	end
end	
