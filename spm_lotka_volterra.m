function [f] = spm_lotka_volterra(x,v,P)
% euations of motion for Lotka Volterra dynamics
% FORMAT spm_lotka_volterra
% x - hidden states
% v - exogenous inputs
% P - 
% 
% returns f = dx/dt = C*S(x) - x;
%              S(x) = 1./(1 + exp(-P*x))
%
% where C determines the order of unstable fixed points visitied in the
% stable heteroclinic channel.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_lotka_volterra.m 2905 2009-03-20 13:00:15Z karl $

% connectivity
%--------------------------------------------------------------------------
x      = spm_vec(x) - 4;
v      = spm_vec(v);
n      = length(x);
m      = length(v);
C      = spm_speye(n,n,0)*(1 - 8) + spm_speye(n,n,1)*(1/2 - 8) + 8;
C(n,1) = 1/2;

% flow
%--------------------------------------------------------------------------
f      = -C*(1./(1 + exp(-P*x))) - x/8;
f      = 2*f;
f(1:m) = f(1:m) + v;