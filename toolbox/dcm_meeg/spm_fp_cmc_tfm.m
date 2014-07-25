function [f] = spm_fp_cmc_tfm(x,u,P,M)
% parameter equations for a neural mass model (canonical microcircuit)
% FORMAT [f] = spm_fp_cmc_tfm(x,u,P,M)
%
% x      - state vector
%   x(:,1) - voltage     (spiny stellate cells)
%   x(:,2) - conductance (spiny stellate cells)
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - voltage     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%   x(:,7) - voltage     (deep pyramidal cells)
%   x(:,8) - conductance (deep pyramidal cells)
%
% f        - dp(t)/dt  = f(x(t),u(t),P,M)
%
% Prior fixed parameter scaling
%
% G  = intrinsic rates
% D  = propagation delays (intrinsic, extrinsic)
% T  = synaptic time constants
% R  = slope of sigmoid activation function
%
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fp_cmc_tfm.m 6123 2014-07-25 17:10:51Z karl $

% Neuronal states (deviations from baseline firing)
%--------------------------------------------------------------------------
%   x(:,1) - voltage     (spiny stellate cells)
%   x(:,2) - conductance (spiny stellate cells)
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - voltage     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%   x(:,7) - voltage     (deep pyramidal cells)
%   x(:,8) - conductance (deep pyramidal cells)
%--------------------------------------------------------------------------
persistent iG nP
if isempty(iG)
    iG  = spm_fieldindices(P,'G');
    nP  = spm_length(P);
end
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
f  = zeros(nP,1);                          % flow
dG = zeros(size(P.G));                     % change in parameters
x  = spm_unvec(x,M.x);                     % neuronal states


% neuronal populations with Voltage-dependent connectivity
%--------------------------------------------------------------------------
i  = [3 1 5];
a  = [32 8 8];                             % potentiation (upper bound)
b  = [2 4 8];                              % decay rate

% NMDA-like Voltage-dependent chnegs i synaptic efficacy
%--------------------------------------------------------------------------
for j = 1:size(P.E,2)
    A       = exp(a(j)*exp(P.E(:,j)).*x(:,i(j))) - 1;
    dG(:,j) = A.*(     exp(P.F(:,j))/2 - P.G(:,j)) - ....
                       b(j)*(P.G(:,j)  - M.Q.G(:,j));
end

f(iG) = dG(:);

