function [y] = spm_gx_lfp(x,u,P,M)
% observer for a neural mass model of erps
% FORMAT [y] = spm_gx_lfp(x,u,P,M)
% x      - state vector
%   x(:,1)  - voltage (spiny stellate cells)
%   x(:,2)  - voltage (pyramidal cells)         +ve
%   x(:,3)  - voltage (pyramidal cells)         -ve
%   x(:,4)  - current (spiny stellate cells)    +ve 
%   x(:,5)  - current (pyramidal cells)         +ve
%   x(:,6)  - current (pyramidal cells)         -ve
%   x(:,7)  - voltage (inhibitory interneurons) +ve
%   x(:,8)  - current (inhibitory interneurons) +ve
%   x(:,9)  - voltage (pyramidal cells)
%   x(:,10) - voltage (inhibitory interneurons) -ve
%   x(:,11) - current (inhibitory interneurons) -ve
%   x(:,12) - voltage (inhibitory interneurons)
%
%   x(:,13) - slow potassium conductance
%
% y        - measured voltage
%__________________________________________________________________________
%
% This is a simplified version of spm_gx_erp
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_lfp.m 1174 2008-02-27 20:22:30Z karl $


% configure state variables
%--------------------------------------------------------------------------
x  = spm_unvec(x,M.x);

% parameterised lead field
%--------------------------------------------------------------------------
L  = spm_erp_L(P,M);
y  = L*x*P.J;

