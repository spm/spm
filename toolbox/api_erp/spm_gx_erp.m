function [y] = spm_gx_erp(x,u,P,M)
% observer for a neural mass model of erps
% FORMAT [y] = spm_gx_erp(x,u,P,M)
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
% y        - measured voltage
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_erp.m 1183 2008-03-03 18:26:05Z karl $

% configure state variables
%--------------------------------------------------------------------------
x  = spm_unvec(x,M.x);

% parameterised lead field
%--------------------------------------------------------------------------
L  = spm_erp_L(P,M);
y  = L*x*P.J';

