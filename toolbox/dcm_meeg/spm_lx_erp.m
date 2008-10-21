function [L] = spm_lx_erp(P,M)
% observer matrix for a neural mass model of erps: y = G*x
% FORMAT [G] = spm_lx_erp(P,M)
% x      - state vector
% G      - where y = G*x; G = dy/dx
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lx_erp.m 2374 2008-10-21 18:52:29Z karl $

% parameterised lead field times source contribution to ECD
%--------------------------------------------------------------------------
L       = spm_erp_L(P,M);                    % lead field per source
L       = kron(P.J,L);                       % lead-field per state
