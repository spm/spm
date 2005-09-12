function [R] = spm_cov2corr(C)
% returns the correlation matrix given the covariance matrix
% FORMAT [R] = spm_cov2corr(C);
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_cov2corr.m 226 2005-09-12 15:00:59Z karl $


%--------------------------------------------------------------------------
n    = length(C);
i    = 1:n;
D    = sparse(i,i,sqrt(1./(diag(C) + eps)));
R    = real(D*C*D);
