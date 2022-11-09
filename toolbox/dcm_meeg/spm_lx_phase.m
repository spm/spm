function [G] = spm_lx_phase(P,M)
% Observation function for phase-coupled oscillators
% FORMAT [G] = spm_lx_phase(P,M)
%
% G     Observations y = Gx
%__________________________________________________________________________
 
% Will Penny
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


Nr = length(P.L);
G  = eye(Nr);
