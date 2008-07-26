function [y] = spm_fx_lz(x,u,P)
% flow for Lorenz attractor
% FORMAT [y] = spm_fx_lz(x,u,P)
% x - state
% u - input
% P - parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_lz.m 1961 2008-07-26 09:38:46Z karl $


% flow for Lorenz attractor
%--------------------------------------------------------------------------
try
    P(3) = P(3)*(1 + u);
end
J    = [-P(1) P(1) 0; (P(3) - x(3)) -1 -x(1); x(2) x(1) P(2)];
y    = J*x;
