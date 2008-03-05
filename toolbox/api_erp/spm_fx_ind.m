function [f,J] = spm_fx_ind(x,u,P,M)
% state equations for a neural mass model of erps
% FORMAT [f,J] = spm_fx_erp(x,u,P,M)
%   x(i,j) - power in the i-th region and j-th frequency mode
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_ind.m 1186 2008-03-05 12:52:57Z karl $


% State: f(x)
%==========================================================================

% Global scaling
%--------------------------------------------------------------------------
K  = 16*exp(P.K);

% neuronal input
%--------------------------------------------------------------------------
U  = 32*P.C*u(:);

% Bilinear form
%--------------------------------------------------------------------------
f  = K*(P.A*x(:) + U);

if nargout == 1, return, end

% Jacobian: J = df(x)/dx
%==========================================================================

% augment with time
%--------------------------------------------------------------------------
J  = K*P.A;





