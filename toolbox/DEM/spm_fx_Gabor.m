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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_Gabor.m 1143 2008-02-07 19:33:33Z spm $

% state variables
%--------------------------------------------------------------------------
f  = u - x; 

