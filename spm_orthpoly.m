function C = spm_orthpoly(N,K)
% Create orthonormal polynomial basis functions
% FORMAT C = spm_orthpoly(N,[K])
% N - dimension
% K - order
%__________________________________________________________________________
% spm_orthpoly creates a matrix for the first few basis functions of an
% orthogonal polynomial expansion
%__________________________________________________________________________
% Copyright (C) 2007 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_orthpoly.m 5900 2014-02-27 21:54:51Z karl $

%--------------------------------------------------------------------------
if nargin == 1, K = N; end
C     = zeros(N,K + 1);
x     = (1:N)';
for i = 0:K
    C(:,i + 1) = x.^i;
end
C = spm_orth(C,'norm');


