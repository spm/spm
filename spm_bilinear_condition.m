function M0 = spm_bilinear_condition(M0)
% conditions a bilinear operator by suppressing positive eigenmodes
% FORMAT M0 = spm_bilinear_condition(M0)
% M0 - bilinear operator
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_bilinear_condition.m 6937 2016-11-20 12:30:40Z karl $
 
% conditions a bilinear operator by suppressing positive eigenmodes
%==========================================================================
 
% remove unstable modes from Jacobian
%--------------------------------------------------------------------------
dfdx  = M0(2:end,2:end);
[u,s] = eig(full(dfdx),'nobalance');
s     = diag(s);
s     = 1j*imag(s) + real(s) - exp(real(s));
 
% replace in bilinear operator
%--------------------------------------------------------------------------
M0(2:end,2:end) = real(u*diag(s)*spm_pinv(u));
