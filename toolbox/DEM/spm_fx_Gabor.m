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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_fx_Gabor.m 1135 2008-02-06 14:15:16Z karl $

% state variables
%--------------------------------------------------------------------------
f  = u - x; 

