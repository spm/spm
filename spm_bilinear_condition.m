function M0 = spm_bilinear_condition(M0,t)
% conditions a bilinear operator by suppressing positive eigenmodes
% FORMAT M0 = spm_bilinear_condition(M0,t)
% M0 - bilinear operator
% t  - time constant of largest eigenmode (secs {default: 32]
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_bilinear_condition.m 6856 2016-08-10 17:55:05Z karl $
 
% conditions a bilinear operator by suppressing positive eigenmodes
%==========================================================================
 
 
% upper bound on unstable modes (32 seconds)
%--------------------------------------------------------------------------
if nargin == 1, t = 16; end
 
% remove unstable modes from Jacobian
%--------------------------------------------------------------------------
dfdx  = M0(2:end,2:end);
[u,s] = eig(full(dfdx),'nobalance');
s     = diag(s);
s     = 1j*imag(s) + real(s) - exp(real(s));
 
% replace in bilinear operator
%--------------------------------------------------------------------------
M0(2:end,2:end) = real(u*diag(s)*spm_pinv(u));
