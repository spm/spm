function C = spm_cov2corr(C)
% Correlation matrix given the covariance matrix
% FORMAT R = spm_cov2corr(C)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


n = length(C);
D = sparse(1:n,1:n,sqrt(1./(diag(C) + eps)));
C = real(D*C);
C = real(C*D);
