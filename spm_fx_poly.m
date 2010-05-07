function [f] = spm_fx_poly(x,v,P)
% Normal (bilinear) form equation of motion
% FORMAT [f] = spm_fx_poly(x,v,P)
% x      - state vector
% v      - exogenous cause
% P      - free parameters 
%
% f      - dx/dt
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_poly.m 3878 2010-05-07 19:53:54Z karl $

% compute Jacobian from blinear terms
%--------------------------------------------------------------------------
x     = spm_vec(x);
J     = P.A;
for i = 1:length(P.B)
    J = J + P.B{i}*x(i);
end
for i = 1:length(P.C)
    J = J + P.C{i}*v(i);
end
f     = J*x;
    
