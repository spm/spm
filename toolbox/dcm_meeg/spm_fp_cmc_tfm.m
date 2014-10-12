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
% f        - dP = h(x(t),u(t),P,M)
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
% $Id: spm_fp_cmc_tfm.m 6234 2014-10-12 09:59:10Z karl $

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
persistent iG nP Ca G
if isempty(iG)
    Ca  = zeros(size(P.G));
    G   = zeros(size(P.G));                    % parameter (deviates)
    iG  = spm_fieldindices(P,'G');
    nP  = spm_length(P);
end
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
f  = zeros(nP,1);                              % flow
x  = spm_unvec(x,M.x);                         % neuronal states
x  = x(:,1:2:end);                             % depolarisation

% neuronal populations with Voltage-dependent connectivity
%==========================================================================
%                  ss sp ii dp                 % neuronal populations
%--------------------------------------------------------------------------
a     = [1 8 2 1]*64;                          % potentiation rate
b     = [4 2 2 1]*4;                           % decay rate

NMDA  = @(x)1./(1 + exp(-x)) - 1/2;            % depolarisation CDF

% NMDA-like Voltage-dependent changes in (recurrent) synaptic efficacy
%--------------------------------------------------------------------------
A     = exp(P.E)*diag(a);
B     = exp(P.F)*diag(b);
dC    = (A.*NMDA(8*x) - Ca).*B;
dG    = Ca.*(2 - G)/2  - G.*B;
Ca    = Ca + dC*M.dt;
G     = G  + dG*M.dt;
f(iG) = G(:);


