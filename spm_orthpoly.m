function C = spm_orthpoly(N,K)
% Create orthonormal polynomial basis functions
% FORMAT C = spm_orthpoly(N,[K])
% N - dimension
% K - order
%__________________________________________________________________________
%
% spm_orthpoly creates a matrix for the first few basis functions of an
% orthogonal polynomial expansion.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2007-2022 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
if nargin == 1, K = N; end
C     = zeros(N,K + 1);
x     = (1:N)';
for i = 0:K
    C(:,i + 1) = x.^i;
end
C = spm_orth(C,'norm');
