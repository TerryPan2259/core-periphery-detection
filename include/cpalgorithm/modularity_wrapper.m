% Modularity wrapper 
%
% @version 0.1.0
% @date 2017/02/19
%
classdef modularity_wrapper < cpabst
	properties
		md;
	end

	methods ( Access = public )
		function obj = modularity_wrapper()
				
			addpath /panfs/panasas01/emat/sk16722/program/lib/comalgorithm
			obj.md = modularity();
			
		end	

		function [C, P, Q, Qs, score, param, cpu_time] = detect( self, G, param )
			[C,Q,Qs,score,param,cpu_time] = self.md.detect( G, param );
			P = sparse(size(C,1),size(C,2),0);
		end
		
		function [Q, Qs, score] = eval( self, G, C, P, param )
			[Q,Qs,score] = self.md.eval(G,C,param);	
		end
			
		function param = initParam(self,param)
			param = self.md.initParam(param);
		end
	end
end
