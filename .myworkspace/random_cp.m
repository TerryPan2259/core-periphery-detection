% Discrete variant of the Borgatti-Everett algorithm 
%
% @version 0.1.0
% @date 2017/02/19
%
% @author Sadamori Kojaku 
classdef random_cp < cpabst

	methods ( Access = public )

		function [C, P, Q, Qs, score, param, cpu_time] = detect( self, G, param )
			% --------------------------------
			% Initialise
			% --------------------------------
			C= []; P = []; Q = -Inf; score = [];
			ts = cputime;
			
			% --------------------------------
			% Generate random partition 
			% --------------------------------
			C = sparse(G.numNode(),1,0); 
			C(rand(G.numNode(),1) < 0.5) = 1;P = 1-C;	
			
			[Q,Qs,score] = self.eval( G, C );	
			cpu_time = cputime-ts;
		end
		
		function [q,qs,score] = eval(self,G,x,varargin)
			q = 0;
			qs = 0;
			score = sparse(G.numNode(),1);
		end
			
		function param = initParam(self,param)
			if ~isfield(param,'numRun');param.numRun = 1;end
			if ~isfield(param,'name') param.name = 'random';end
			if ~isfield(param,'disp') param.disp = false;end
		end
	end
	
end
