function [f] = spm_mc_fxa_5(x,v,a,P)
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
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% physical flow
%--------------------------------------------------------------------------
f   = x;
dx  = spm_fx_mountaincar([x.x; x.v],v,a,P)/2;
f.x = dx(1);
f.v = dx(2);

% physiological flow
%--------------------------------------------------------------------------
f.q = exp(-(x.x - 1)^2*32) - x.q/32;
f.p = (x.q - x.p)/64;
