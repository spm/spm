function p = spm_ness_Sp2p(Sp,x,K)
% converts a density into polynomial potential parameters  
% FORMAT p = spm_ness_Sp2p(Sp,x,[K])
%--------------------------------------------------------------------------
% Sp   - Polynomial coefficients or parameters of log density
% x{i} - support (sample points): i = 1,...,N
% K    - Order of polynomial expansion (K = 3 corresponds to quadratic)
%
% p    - probability density
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_cond.m 8097 2021-04-24 20:28:27Z karl $

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

return

