function [f,J] = spm_fx_ind(x,u,P,M)
% State equations for a neural mass model of erps
% FORMAT [f,J] = spm_fx_erp(x,u,P,M)
%   x(i,j) - power in the i-th region and j-th frequency mode
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% Jacobian: J = df(x)/dx
%==========================================================================
J  = P.A;

% State: f(x)
%==========================================================================
f  = J*x(:) + P.C*u(:);
