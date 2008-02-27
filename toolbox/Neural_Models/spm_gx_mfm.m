function [y] = spm_gx_lfp(x,u,P,M)
% observer for a neural mass model of erps
% FORMAT [y] = spm_gx_lfp(x,u,P,M)
%
% x      - state vector
% x{i,j} - states of i-th source and j-th population
% 
% x{i,1} - excitatory spiny stellate cells (input cells)
% x{i,2} - inhibitory interneurons
% x{i,3} - excitatory pyramidal cells      (output cells)
%
% x{i,j}.V  - voltage
% x{i,j}.gE - conductance (excitatory)
% x{i,j}.gI - conductance (inhibitory)
%
% y        - measured voltage
%__________________________________________________________________________
%
% This is a simplified version of spm_gx_erp
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%


% get dimensions and configure state variables
%--------------------------------------------------------------------------
n  = length(P.A{1});
x  = spm_vec(x{:,end}.V);

% parameterised lead field
%--------------------------------------------------------------------------
L  = spm_erp_L(P,M);
y  = L*x;


