function [f,J] = spm_fx_ind(x,u,P,M)
% state equations for a neural mass model of erps
% FORMAT [f,J] = spm_fx_erp(x,u,P,M)
%   x(i,j) - power in the i-th region and j-th frequency mode
%             NB: the first state is actually time but this hidden here.
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_ind.m 1040 2007-12-21 20:28:30Z karl $

% check input u = f(t,P,M)
%--------------------------------------------------------------------------
try
    fu  = M.fu;
catch
    fu  = 'spm_ind_u';
end

% get dimensions and configure state variables
%--------------------------------------------------------------------------
t  = x(1);                       % peristimulus time (sec)
x  = x(2:end);                   % neuronal states

% dfdx = [] if t exceeds trial duration (invoking a return to initial state)
%--------------------------------------------------------------------------
if nargout == 1 & (t - M.dur) > 1e-6, f = []; return, end

% effective extrinsic connectivity
%--------------------------------------------------------------------------
for i = 1:M.m
	  P.A = P.A + u(i)*P.B{i};
end

% Global scaling
%--------------------------------------------------------------------------
K  = 16*exp(P.K);

% input
%--------------------------------------------------------------------------
U  = feval(fu,t,P,M);
U  = 32*P.C*U;

% State: f(x)
%==========================================================================

% Bilinear form
%--------------------------------------------------------------------------
f  = K*(P.A*x + U);

% augment with time
%--------------------------------------------------------------------------
f  = [1; f(:)];

if nargout == 1, return, end

% Jacobian: J = df(x)/dx
%==========================================================================

% augment with time
%--------------------------------------------------------------------------
J  = blkdiag(0,K*P.A);





