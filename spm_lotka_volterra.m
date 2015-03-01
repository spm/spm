function [f] = spm_lotka_volterra(x,v,P)
% equations of motion for Lotka-Volterra dynamics
% FORMAT [f] = spm_lotka_volterra(x,v,P)
% FORMAT [f] = spm_lotka_volterra(x,v)
% FORMAT [P] = spm_lotka_volterra(n)
%
% x     - hidden states
% v     - parameter of P.f
% P     - lateral connectivity
%
% returns f = dx/dt = P*S(x) - x/8 + 1;
%              S(x) = 1./(1 + exp(-x))
%
% where P determines the order of unstable fixed points visited in the
% stable heteroclinic channel. If P is not specified it will be computed
% using v. If x is a scalar a matrix of size x (P) is returned (with v = 1).
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_lotka_volterra.m 6353 2015-03-01 11:52:49Z karl $


% intialise
%==========================================================================
try, v; catch, v = 1; end

% just return connectivity
%--------------------------------------------------------------------------
if isscalar(x) && isnumeric(x)
    n  = x;
    P  = spm_speye(n,n,-1) - spm_speye(n,n,1); P(n,1) = -1; P(1,n) = 1;
    f  = v*P + speye(n,n) - 1;
    return
end

% check for parameters of succession
%--------------------------------------------------------------------------
if nargin < 3
    P  = spm_lotka_volterra(spm_length(x),v);
end


% flow
%==========================================================================
f   = P*(1./(1 + exp(-x))) - x/8 + 1;
    

