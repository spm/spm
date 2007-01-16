function [C] = spm_cov2corr(C)
% returns the correlation matrix given the covariance matrix
% FORMAT [R] = spm_cov2corr(C);
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_cov2corr.m 716 2007-01-16 21:13:50Z karl $


%--------------------------------------------------------------------------
n    = length(C);
D    = sparse(1:n,1:n,sqrt(1./(diag(C) + eps)));
C    = real(D*C);
C    = real(C*D);
