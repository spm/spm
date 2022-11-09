function p = spm_ness_Sp2p(Sp,x,K)
% Convert a density into polynomial potential parameters  
% FORMAT p = spm_ness_Sp2p(Sp,x,[K])
%--------------------------------------------------------------------------
% Sp   - Polynomial coefficients or parameters of log density
% x{i} - support (sample points): i = 1,...,N
% K    - Order of polynomial expansion (K = 3 corresponds to quadratic)
%
% p    - probability density
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% defaults
%--------------------------------------------------------------------------
if nargin < 3, K = 3; end

for i = 1:numel(x)
    N(i) = numel(x{i});
end

% ensure p is a proper density
%--------------------------------------------------------------------------
p = spm_softmax(spm_polymtx(x,K)*Sp);
p = reshape(p,N);
