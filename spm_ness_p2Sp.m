function Sp = spm_ness_p2Sp(p,x,K)
% Convert a density into polynomial potential parameters  
% FORMAT Sp = spm_ness_p2Sp(p,x,K))
%--------------------------------------------------------------------------
% p    - probability density
% x{i} - support (sample points): i = 1,...,N
% K    - Order of polynomial expansion (K = 3 corresponds to quadratic)
%
% Sp   - Polynomial coefficients or parameters of log density
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% defaults
%--------------------------------------------------------------------------
if nargin < 3, K = 3; end

% ensure p is a proper density
%--------------------------------------------------------------------------
p       = p(:) + eps;
p       = p/sum(p);

% and corresponding polynomial coefficients
%--------------------------------------------------------------------------
b       = spm_polymtx(x,K);
Sp      = b\log(p);
