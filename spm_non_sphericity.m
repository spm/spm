function [xVi] = spm_non_sphericity(xVi)
% return error covariance constraints for basic ANOVA designs
% FORMAT [xVi] = spm_non_sphericity(xVi)
%
% required fields:
% xVi.I    - n x 4 matrix of factor level indicators
%              I(n,i) is the level of factor i for observation n
% xVi.var  - 1 x 4 vector of flags
%              var(i) = 1; levels of factor i have unequal variances
% xVi.dep  - 1 x 4 vector of flags
%              dep(i) = 1; levels of factor i are dependent
%
% Output:
% xVi.Vi   -  cell of covariance components
%
% See also; spm_Ce.m & spm_spm_ui.m
%___________________________________________________________________________
% %W% Karl Friston %E%


% Identify replication factor (r)
%---------------------------------------------------------------------------
I     = xVi.I;
r     = find((max(I) > 1) & ~xVi.var & ~xVi.dep);
try
	r = r(1);
catch
	error('There are no replications to estimate covariances')
end

% create covariance components Q{:}
%===========================================================================
[n f] = size(I);			% # observations, % # Factors
R     = sparse([1:n],I(:,r),1);		% main effect of replication factor
RR    = R*R';
Q     = {};

% unless all repeated measures are identically distributed
%---------------------------------------------------------------------------
if ~any(xVi.var)
	Q{end + 1} = speye(n,n);
end

for i = 1:f

	% unequal variances of repeated measure levels
	%-------------------------------------------------------------------
	nL    = max(I(:,i));
	if xVi.var(i)

		% add variance component for level j of factor i
		% i.e. variance due to Fi(j) x R
		%-----------------------------------------------------------
		for j = 1:nL
			u  = I(:,i) == j;
			q  = spdiags(u,0,n,n);
			Q{end + 1} = q*RR*q;
		end
	end

	% dependencies among repeated measure levels
	%-------------------------------------------------------------------
	if xVi.dep(i)

		% add covariance component for levels j & k of factor i
		% i.e. covariance between Fi(j) x R & Fi(k) x R
		%-----------------------------------------------------------
		for j = 1:nL
		for k = (j + 1):nL
			u  = I(:,i) == j;
			v  = I(:,i) == k;
			p  = spdiags(u,0,n,n);
			q  = spdiags(v,0,n,n);
			Q{end + 1} = p*RR*q + q*RR*p;
		end
		end
	end
end

% set Q in non-sphericity structure
%---------------------------------------------------------------------------
xVi.Vi = Q;
