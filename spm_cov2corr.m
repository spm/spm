function [R] = spm_cov2corr(C)
% returns the correlation matrix given the covariance matrix
% FORMAT [R] = spm_cov2corr(C);
%___________________________________________________________________________
% @(#)spm_cov2corr.m	2.1 Karl Friston 05/01/10

%---------------------------------------------------------------------------
n    = length(C);
i    = 1:n;
D    = sparse(i,i,sqrt(1./(diag(C) + eps)));
R    = D*C*D;
