function [f] = mci_nmm_fx_delay (x,u,P,M)
% State equations for a neural mass model of erps with first order delays
% FORMAT [f] = mci_nmm_fx_delay (x,u,P,M)
%
% x      - state vector
%   x(:,1) - voltage (spiny stellate cells)
%   x(:,2) - voltage (pyramidal cells) +ve
%   x(:,3) - voltage (pyramidal cells) -ve
%   x(:,4) - current (spiny stellate cells)    depolarizing
%   x(:,5) - current (pyramidal cells)         depolarizing
%   x(:,6) - current (pyramidal cells)         hyperpolarizing
%   x(:,7) - voltage (inhibitory interneurons)
%   x(:,8) - current (inhibitory interneurons) depolarizing
%   x(:,9) - voltage (pyramidal cells)
%
% f        - dx(t)/dt  = f(x(t))
%
% Prior fixed parameter scaling [Defaults]
%
%  M.pF.E = [32 16 4];           % extrinsic rates (forward, backward, lateral)
%  M.pF.H = [1 4/5 1/4 1/4]*128; % intrinsic rates (g1, g2 g3, g4)
%  M.pF.D = [2 16];              % propogation delays (intrinsic, extrinsic)
%  M.pF.G = [4 32];              % receptor densities (excitatory, inhibitory)
%  M.pF.T = [8 16];              % synaptic constants (excitatory, inhibitory)
%  M.pF.R = [1 1/2];             % parameter of static nonlinearity
%
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_nmm_fx_delay.m 6548 2015-09-11 12:39:47Z will $

if ~isstruct(P)
%     if isfield(M,'V')
%         % Project params out of subspace
%         P=M.V*P(:);
%     end
    P=spm_unvec(P,M.pE);
end

% get dimensions and configure state variables
%--------------------------------------------------------------------------
n = length(P.A{1});         % number of sources
x = spm_unvec(x,M.x);       % neuronal states

% [default] fixed parameters
%--------------------------------------------------------------------------
E = [1 1/2 1/8]*32;         % extrinsic rates (forward, backward, lateral)
G = [1 4/5 1/4 1/4]*128;    % intrinsic rates (g1 g2 g3 g4)
D = [2 16];                 % propagation delays (intrinsic, extrinsic)
H = [4 32];                 % receptor densities (excitatory, inhibitory)
T = [8 16];                 % synaptic constants (excitatory, inhibitory)
R = [2 1]/3;                % parameters of static nonlinearity

% [specified] fixed parameters
%--------------------------------------------------------------------------
if isfield(M,'pF')
    try, E = M.pF.E; end
    try, G = M.pF.H; end
    try, D = M.pF.D; end
    try, H = M.pF.G; end
    try, T = M.pF.T; end
    try, R = M.pF.R; end
end

% test for free parameters on intrinsic connections
%--------------------------------------------------------------------------
try
    G = G.*exp(P.H);
end
G     = ones(n,1)*G;

% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1}  = exp(P.A{1})*E(1);
A{2}  = exp(P.A{2})*E(2);
A{3}  = exp(P.A{3})*E(3);
C     = exp(P.C);

% intrinsic connectivity and parameters
%--------------------------------------------------------------------------
Te    = T(1)/1000*exp(P.T(:,1));         % excitatory time constants
Ti    = T(2)/1000*exp(P.T(:,2));         % inhibitory time constants
He    = H(1)*exp(P.G(:,1));              % excitatory receptor density
Hi    = H(2)*exp(P.G(:,2));              % inhibitory receptor density

% delays
%--------------------------------------------------------------------------
De = D(2).*exp(P.D)/1000;                % extrinsic
Di = D(1)/1000;                          % intrinsic

% presynaptic i/p from delayed voltages using first order Taylor approx
% v(t-tau) = v(t) - tau dv/dt
%--------------------------------------------------------------------------
pyr_int = presynaptic (x(:,9)-Di*(x(:,5)-x(:,6)),P,R);  % Intrinsic delays
stl_int = presynaptic (x(:,1)-Di*x(:,4),P,R);
inh_int = presynaptic (x(:,7)-Di*x(:,8),P,R);
pyr_ext = presynaptic (x(:,9)-De*(x(:,5)-x(:,6)),P,R);  % Extrinsic delays

% input
%==========================================================================
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U = u(:)*64;
    
else
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:)*2;
end

% State: f(x)
%==========================================================================

%f=zeros(size(M.x));

% Supragranular layer (inhibitory interneurons): Voltage & depolarizing current
%--------------------------------------------------------------------------
f(:,7) = x(:,8);
f(:,8) = (He.*((A{2} + A{3})*pyr_ext + G(:,3).*pyr_int) - 2*x(:,8) - x(:,7)./Te)./Te;

% Granular layer (spiny stellate cells): Voltage & depolarizing current
%--------------------------------------------------------------------------
f(:,1) = x(:,4);
f(:,4) = (He.*((A{1} + A{3})*pyr_ext + G(:,1).*pyr_int + U) - 2*x(:,4) - x(:,1)./Te)./Te;

% Infra-granular layer (pyramidal cells): depolarizing current
%--------------------------------------------------------------------------
f(:,2) = x(:,5);
f(:,5) = (He.*((A{2} + A{3})*pyr_ext + G(:,2).*stl_int) - 2*x(:,5) - x(:,2)./Te)./Te;

% Infra-granular layer (pyramidal cells): hyperpolarizing current
%--------------------------------------------------------------------------
f(:,3) = x(:,6);
f(:,6) = (Hi.*G(:,4).*inh_int - 2*x(:,6) - x(:,3)./Ti)./Ti;

% Infra-granular layer (pyramidal cells): Voltage
%--------------------------------------------------------------------------
f(:,9) = x(:,5) - x(:,6);

%f      = spm_vec(f);
f=f(:);

end

function [S] = presynaptic (x,P,R)

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R     = R.*exp(P.S);
S     = 1./(1 + exp(-R(1)*(x - R(2)))) - 1./(1 + exp(R(1)*R(2)));

end






