function [f] = spm_SHC_fx(x,v,P,varargin)
% equations of motion for Lotka-Volterra dynamics
% FORMAT [f] = spm_SHC_fx(x,v,P)
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
% $Id: spm_SHC_fx.m 3265 2009-07-10 14:02:22Z karl $



% intialise
%==========================================================================

% SHC states 
%--------------------------------------------------------------------------
f.x  = P.f*spm_phi(x.x) - x.x/32 + 1;
f.c  = (spm_phi(x.x) - x.c)/64;

f    = spm_vec(f)*16;




