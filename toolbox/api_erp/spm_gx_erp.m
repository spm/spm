function [y] = spm_gx_erp(x,u,P)
% observer for a neural mass model of erps
% FORMAT [y] = spm_gx_erp(x,u,P)
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
%___________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%___________________________________________________________________________
% %W% Karl Friston %E%

% get dimensions and configure state variables
%---------------------------------------------------------------------------
x     = x(2:end);
m     = length(u);  				% number of inputs
n     = length(x)/9;  				% number of regions
x     = reshape(x,n,9);             % x

% configure parameters: Af,...
%---------------------------------------------------------------------------
P    = spm_erp_pack(P,m,n);

% output (sources if Lead field is empty)
%===========================================================================
if ~length(P.L)
    y = x(:,9);
else
    y = P.L*(x(:,9).*P.K);
end
