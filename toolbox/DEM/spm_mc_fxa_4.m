function [f] = spm_mc_fxa_4(x,v,a,P)
% equations of motion for the mountain car problem
% problem
% FORMAT [f] = spm_mc_fxa_4(x,v,a,P)
%
% x   - hidden states
% v   - exogenous inputs
% a   - action
% P   - parameters for mountain car
%
% returns f = dx/dt 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mc_fxa_4.m 3333 2009-08-25 16:12:44Z karl $
 
 
% physical flow
%--------------------------------------------------------------------------
dx  = spm_fx_mountaincar([x.x; x.v],v,a,P)/2;
f.x = dx(1);
f.v = dx(2);

% physiological flow
%--------------------------------------------------------------------------
A   = exp(-(x.x - 1).^2*32);
f.p = A - x.p/32;
f.d = (x.p - x.d)/64;
