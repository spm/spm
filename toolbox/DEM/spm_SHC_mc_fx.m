function [f] = spm_SHC_mc_fx(x,v,P)
% equations of motion for SHC dynamics of the mountian car
% problem
% FORMAT [f] = spm_SHC_mc_fx(x,v,P)
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
% $Id: spm_SHC_mc_fx.m 3113 2009-05-11 15:25:13Z karl $



% intialise
%==========================================================================
try, P.k; catch, P.k = 1; end
try, P.l; catch, P.l = 1; end


% SHC states and flow
%--------------------------------------------------------------------------
f.x  = spm_lotka_volterra(x.x,x.v,P);

% and associated attractors
%--------------------------------------------------------------------------
x.e  = exp(4*x.x);
P.g  = P.g*x.e/sum(x.e);
f.v  = [ x.v(2);
         P.g(3)*(P.g(1) - x.v(1)) - x.v(2)*P.g(2)];



