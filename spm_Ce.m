function [C] = spm_Ce(v,h)
% return error covariance contraints for serially correlated data
% FORMAT [C] = spm_Ce(v,h)
% v  - (1 x l) v(i) = number of obervations for ith block
% h  - length of expansion = length(C) (default h = 1)
%
% C  - {i x lh} structure of block diagonal (lv x lv)matrices
%	C{1:l + 0*l} speye(v,v) - i.i.d. errors [coeficient = 1/(0*e)]
%	C{1:l + 1*l} AR(1) with coefficient 1/(1*e)
%	C{1:l + 2*l} AR(1) with coefficient 1/(2*e)
%
%		.....
%
%	C{1:l + (h - 1)}
%___________________________________________________________________________

% defaults
%---------------------------------------------------------------------------
if nargin == 1
	h = 1;
end
l      = length(v);

% create constraints on Cov{e} - C1
%===========================================================================
Q      = [];
q      = 6;
for  i = 1:h
	k      = exp(-[0:q]'/(i - 1 + eps));
	k      = conv(k,flipud(k));
	if length(Q)
		k = k - Q*(pinv(Q)*k);
	end
	Q(:,i) = k;
end


% C
%---------------------------------------------------------------------------
C      = {};
for  i = 1:h
	for j = 1:l
		k           = [1:v(j)] + sum(v(1:(j - 1)));
		K           = spdiags(ones(v(j),1)*Q(:,i)',[-q:q],v(j),v(j));
		C{end + 1}  = sparse(sum(v),sum(v));
		C{end}(k,k) = sparse(K.*(abs(K) > 1e-2));
	end
end

