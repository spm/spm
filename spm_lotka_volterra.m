function [f] = spm_lotka_volterra(x,v,P)
% equations of motion for Lotka-Volterra dynamics
% FORMAT [f] = spm_lotka_volterra(x,v,P)
%
% x   - hidden states
% v   - exogenous inputs
% P.f - lateral connectivity
% P.k - rate [default 1]
%
% returns f = dx/dt = P.f*S(x) - x/8 + 1;
%              S(x) = 1./(1 + exp(-x))
%
% where C determines the order of unstable fixed points visited in the
% stable heteroclinic channel.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_lotka_volterra.m 3113 2009-05-11 15:25:13Z karl $



% intialise
%==========================================================================
try, P.k; catch, P.k = 1; end
try, P.l; catch, P.l = 1; end

try
    
    % SHC states 
    %----------------------------------------------------------------------
    f  = P.f*(1./(1 + exp(-x))) - x/8 + 1;
    f  = P.k*f;

catch

    % SHC states and flow to point attractors in P.g
    %----------------------------------------------------------------------
    x.e  = exp(x.x);
    f.x  = P.f*(1./(1 + exp(-x.x))) - x.x/8 + 1;
    f.v  = P.g*x.e - x.v*sum(x.e);

    f.x  = P.k*f.x;
    f.v  = P.l*f.v;

end

