function [R] = spm_cov2corr(C)
% returns the correlation matrix given the covariance matrix
% FORMAT [R] = spm_cov2corr(C);
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$


%---------------------------------------------------------------------------
n    = length(C);
i    = 1:n;
D    = sparse(i,i,sqrt(1./(diag(C) + eps)));
R    = D*C*D;
