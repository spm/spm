function [C] = spm_cov2corr(C)
% returns the correlation matrix given the covariance matrix
% FORMAT [R] = spm_cov2corr(C);
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cov2corr.m 1143 2008-02-07 19:33:33Z spm $


%--------------------------------------------------------------------------
n    = length(C);
D    = sparse(1:n,1:n,sqrt(1./(diag(C) + eps)));
C    = real(D*C);
C    = real(C*D);
