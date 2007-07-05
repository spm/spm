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
% %W% Karl Friston %E%


% get dimensions and configure state variables
%--------------------------------------------------------------------------
n  = length(P.A{1});
x  = reshape(x,n,13);

% parameterised lead field ECD (pre-multiplied by projector M.E)
%--------------------------------------------------------------------------
if isfield(P,'Lpos')
    L = spm_erp_L(P,M);
else
    L = P.L;
end
try
    y = L*x*P.M;
catch
    y = L*x(:,9);
end

