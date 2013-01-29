function M0 = spm_bilinear_condition(M0,N,dt)
% conditions a bilinear operator by suppressing positive eigenmodes
% FORMAT M0 = spm_bilinear_condition(M0)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_bilinear_condition.m 5219 2013-01-29 17:07:07Z spm $

% conditions a bilinear operator by suppressing positive eigenmodes
%==========================================================================


% regulariser (1/8 of kernel support)
%--------------------------------------------------------------------------
if nargin == 1
    t = 16;
else
    t = 16/(N*dt);
end

% remove unstable modes from Jacobian
%--------------------------------------------------------------------------
dfdx  = M0(2:end,2:end);
[u,s] = eig(full(dfdx));
S     = diag(s);
i     = find(real(S) > -t);
S(i)  = sqrt(-1)*imag(S(i)) - log(exp(t) + exp(real(-S(i))));

% replace in bilinear operator
%--------------------------------------------------------------------------
M0(2:end,2:end) = real(u*diag(S)*spm_pinv(u));
