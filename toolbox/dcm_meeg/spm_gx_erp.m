function [y] = spm_gx_erp(x,~,P,M)
% observer for a neural mass model of event related potentials
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
% $Id: spm_gx_erp.m 4814 2012-07-30 19:56:05Z karl $

% parameterised lead field times [perturbations] of states
%--------------------------------------------------------------------------
L  = spm_lx_erp(P,M.dipfit);
 
% project onto spatial modes
%--------------------------------------------------------------------------
if isfield(M,'U')
    L = M.U'*L;
end
 
% sensor-space response
%--------------------------------------------------------------------------
y  = L*(spm_vec(x) - spm_vec(M.x));

