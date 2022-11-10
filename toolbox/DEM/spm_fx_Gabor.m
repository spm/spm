function [f] = spm_fx_Gabor(x,u,P)
% state equation for Gabor patches
% FORMAT [f] = spm_fx_Gabor(x,u,P)
% x      - state vector
%   x(1) - position
%   x(2) - amplitude
%   x(3) - dispersion
% u      - input
%   u(1) - position   (forcing)
%   u(2) - amplitude  (forcing)
%   u(3) - dispersion (forcing)
% f      - dx/dt
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% state variables
%--------------------------------------------------------------------------
f  = u - x; 

