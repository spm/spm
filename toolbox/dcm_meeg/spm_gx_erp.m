function [y] = spm_gx_erp(x,u,P,M)
% observer for a neural mass model of erps
% FORMAT [y] = spm_gx_erp(x,u,P,M)
% x      - state vector
% y      - measured voltage y = L*x(:)
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_erp.m 2393 2008-10-23 14:58:50Z karl $

% parameterised lead field times [perturbations] of states
%--------------------------------------------------------------------------
L  = spm_lx_erp(P,M);
y  = L*(spm_vec(x) - spm_vec(M.x));

