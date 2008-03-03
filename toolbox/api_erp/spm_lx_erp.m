function [L] = spm_lx_erp(P,M)
% observer matrix for a neural mass model of erps: y = G*x
% FORMAT [G] = spm_lx_erp(P,M)
% x      - state vector
%   x(:,1) - voltage (spiny stellate cells)
%   x(:,2) - voltage (pyramidal cells) +ve
%   x(:,3) - voltage (pyramidal cells) -ve
%   x(:,4) - current (spiny stellate cells)    depolarizing
%   x(:,5) - current (pyramidal cells)         depolarizing
%   x(:,6) - current (pyramidal cells)         hyerpolarizing
%   x(:,7) - voltage (inhibitory interneurons)
%   x(:,8) - current (inhibitory interneurons) depolarizing
%   x(:,9) - voltage (pyramidal cells)
%
% G        - where y = G*x; G = dy/dx
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lx_erp.m 1183 2008-03-03 18:26:05Z karl $

% parameterised lead field times source contribution to ECD
%--------------------------------------------------------------------------
L = kron(P.J,spm_erp_L(P,M));
