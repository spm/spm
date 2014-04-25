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
% $Id: spm_fp_cmc_tfm.m 5966 2014-04-25 14:37:59Z karl $

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
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
f  = spm_unvec(zeros(size(P)),M.pE);  % flow of parameters
x  = spm_unvec(x,M.x);                % neuronal states
P  = spm_unvec(P,M.pE);               % neuronal parameters


% neuronal populations with Voltage-dependent connectivity
%--------------------------------------------------------------------------
A        = exp(32*exp(P.E(:,1)).*x(:,3));     % NMDA
B        = 4*exp(P.F(:,1));                   % ddecay
f.G(:,1) = (A - 1).*(1 - P.G(:,1))  - B.*P.G(:,1);


% f.G(:,1) = exp(32*exp(P.E(:,1)).*x(:,3)) - 1 - 4*exp(P.F(:,1)).*P.G(:,1);
% f.G(:,2) = exp( 8*exp(P.E(:,2)).*x(:,5)) - 1 - 4*exp(P.F(:,2)).*P.G(:,2);

% vectorise
%--------------------------------------------------------------------------
f        = spm_vec(f);


