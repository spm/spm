function [C] = spm_Ce(v,a)
% return error covariance constraints for serially correlated data
% FORMAT [C] = spm_Ce(v,a)
% v  - (1 x l) v(i) = number of obervations for ith block
% a  - AR coeficient expansion point  (default a = [])
% 
%  C{1} = h(1)*AR(a)
%  C{2} = h(1)*AR(a) + h(2)*dAR(a)/da(1);
%  C{3} = h(1)*AR(a) + h(2)*dAR(a)/da(1) + h(3)*dAR(a)/da(2);
%
% See also; spm_Q.m
%___________________________________________________________________________
% %W% Karl Friston %E%


% defaults
%---------------------------------------------------------------------------
if nargin == 1
	a = [];
end


% create blocks
%---------------------------------------------------------------------------
C    = {};
l    = length(v);
n    = sum(v);
k    = 0;
if l > 1
	for i = 1:l
		dCdh  = spm_Ce(v(i),h);
		for j = 1:length(dCdh)
			[x y q]    = find(dCdh{j});
			x          = x    + k;
			y          = y    + k;
			k          = v(i) + k;
			C{end + 1} = sparse(x,y,q,n,n);
		end
	end
else
	% dCdh
	%==================================================================
	C{1}  = spm_Q(a,v);
	dCdh  = spm_diff('spm_Q',a,v,1);
	for i = 1:length(a)
		C{i + 1} = reshape(dCdh(:,i),v,v);
	end

end
