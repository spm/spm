function [L] = spm_lx_dem(P,M)
% Observer matrix for a neural mass model of erps: y = G*x
% FORMAT [G] = spm_lx_dem(P,M)
% x      - state vector
% G      - where y = G*x; G = L*J
%          L = dy/dsource
%          J = dsource/dstate
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% parameterised lead field times source contribution to ECD
%--------------------------------------------------------------------------
L = spm_erp_L(P,M.dipfit)*P.J;                       % lead field per state
