function [xVi] = spm_non_sphericity(xVi)
% return error covariance constraints for basic ANOVA designs
% FORMAT [xVi] = spm_non_sphericity(xVi)
%
% required fields:
% xVi.I    - n x 4 matrix of factor level indicators
%              I(n,i) is the level of factor i for observation n
% xVi.sF   - 1 x 4 cellstr containing factor names 
%              sF{i} is the name of factor i
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


% check there is a replication factor
%---------------------------------------------------------------------------
if ~any(xVi.var | xVi.dep)
	warning('There are no replications to estimate covariances')
end

% create covariance components Q{:}
%===========================================================================

% create variance components
%---------------------------------------------------------------------------
[n f] = size(xVi.I);
Q     = {};
for i = 1:f
	if xVi.var(i)

		% add variance component for this level
		%-----------------------------------------------------------
		nL    = max(xVi.I(:,i));
		for j = 1:nL
			u          = find(xVi.I(:,i) == j);
			q          = sparse(u,u,1,n,n);
			Q{end + 1} = q;
		end
	end
end

% unless all repeated measures are identically distributed
%---------------------------------------------------------------------------
if length(Q) == 0
	Q{1}  = speye(n,n);
end

% independently distributed assumptions
%---------------------------------------------------------------------------
for i = 1:f
	if xVi.dep(i)

		% add covariance component for these levels
		%-----------------------------------------------------------
		nL    = max(xVi.I(:,i));
		for j = 1:nL
		for k = (j + 1):nL
			u          = find(xVi.I(:,i) == j);
			v          = find(xVi.I(:,i) == k);
			q          = sparse(u,v,1,n,n);
			Q{end + 1} = q + q';
		end
		end
	end
end

% set Q in non-sphericity structure
%---------------------------------------------------------------------------
xVi.Vi = Q;
