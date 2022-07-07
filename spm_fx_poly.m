function [f] = spm_fx_poly(x,v,P)
% Normal (bilinear) form equation of motion
% FORMAT [f] = spm_fx_poly(x,v,P)
% x      - state vector
% v      - exogenous cause
% P      - free parameters
%
% f      - dx/dt
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging

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
