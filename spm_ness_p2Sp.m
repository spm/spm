function Sp = spm_ness_p2Sp(p,x,K)
% converts a density into polynomial potential parameters  
% FORMAT Sp = spm_ness_p2Sp(p,x,K))
%--------------------------------------------------------------------------
% p    - probability density
% x{i} - support (sample points): i = 1,...,N
% K    - Order of polynomial expansion (K = 3 corresponds to quadratic)
%
% Sp   - Polynomial coefficients or parameters of log density
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_cond.m 8097 2021-04-24 20:28:27Z karl $

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

return

