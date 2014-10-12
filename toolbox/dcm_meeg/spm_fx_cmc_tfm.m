function [f,J,Q] = spm_fx_cmc_tfm(x,u,P,M)
% state equations for a neural mass model (canonical microcircuit)
% FORMAT [f,J,D] = spm_fx_cmc_tfm(x,u,P,M)
% FORMAT [f,J]   = spm_fx_cmc_tfm(x,u,P,M)
% FORMAT [f]     = spm_fx_cmc_tfm(x,u,P,M)
% x      - state vector
%   x(:,1) - voltage     (spiny stellate cells)
%   x(:,2) - conductance (spiny stellate cells)
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - current     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%   x(:,7) - voltage     (deep pyramidal cells)
%   x(:,8) - conductance (deep pyramidal cells)
% u        - exogenous input
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
% D        - delay operator dx(t)/dt = f(x(t - d))
%                                    = D(d)*f(x(t))
%
% Prior fixed parameter scaling [Defaults]
%
% E  = (forward and backward) extrinsic rates 
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
% $Id: spm_fx_cmc_tfm.m 6234 2014-10-12 09:59:10Z karl $
 
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
x  = spm_unvec(x,M.x);                % neuronal states
n  = size(x,1);                       % number of sources


% [default] fixed parameters
%--------------------------------------------------------------------------
E  = [4 2 2 2]*200;                   % extrinsic (forward and backward)  
T  = [256 128 16 8];                  % synaptic rate constants
R  = 1;                               % gain of activation function
B  = 0;                               % baseline firing

% Extrinsic connections
%--------------------------------------------------------------------------
% ss = spiny stellate
% sp = superficial pyramidal
% dp = deep pyramidal
% ii = inhibitory interneurons
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*E(1);          % forward  connections (sp -> ss)
A{2} = exp(P.A{2})*E(2);          % forward  connections (sp -> dp)
A{3} = exp(P.A{3})*E(3);          % backward connections (dp -> sp)
A{4} = exp(P.A{4})*E(4);          % backward connections (dp -> ii)

% detect and reduce the strength of reciprocal (lateral) connections
%--------------------------------------------------------------------------
for i = 1:length(A)
    L    = (A{i} > exp(-8)) & (A{i}' > exp(-8));
    A{i} = A{i}./(1 + 4*L);
end

% input connections
%--------------------------------------------------------------------------
C    = exp(P.C);
 
% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
F    = 1./(1 + exp(-R*x(:,1:2:end) + B));    % firing rate
S    = F - 1/(1 + exp(B));                   % deviation 

% input
%==========================================================================
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U   = u(:)*256;
    M.m = size(U,1);
    
else
    % exogenous input
    %----------------------------------------------------------------------
    U   = C*u(:)*4;
    M.m = size(C,2);
    
end

 
% time constants and intrinsic connections
%==========================================================================
T      = ones(n,1)*T;
i      = 1:size(P.T,2);
T(:,i) = T(:,i).*exp(P.T);

% extrinsic connections
%--------------------------------------------------------------------------
% forward  (i)   2  sp -> ss (+ve)
% forward  (ii)  1  sp -> dp (+ve)
% backward (i)   2  dp -> sp (-ve)
% backward (ii)  1  dp -> ii (-ve)
%--------------------------------------------------------------------------

% Neuronal states (deviations from baseline firing)
%--------------------------------------------------------------------------
%   S(:,1) - voltage     (spiny stellate cells)
%   S(:,2) - conductance (spiny stellate cells)
%   S(:,3) - voltage     (superficial pyramidal cells)
%   S(:,4) - conductance (superficial pyramidal cells)
%   S(:,5) - current     (inhibitory interneurons)
%   S(:,6) - conductance (inhibitory interneurons)
%   S(:,7) - voltage     (deep pyramidal cells)
%   S(:,8) - conductance (deep pyramidal cells)
%--------------------------------------------------------------------------
%     ss sp ii dp   % intrinsic connections
%--------------------------------------------------------------------------
g  = [-8 -4 -4  0;  % ss
       4 -8 -2  0;  % sp
       4  2 -4  2;  % ii
       0  0 -2 -4]; % dp

g  = g*200*exp(P.S);

% intrinsic connections to be optimised (only the first is modulated)
%--------------------------------------------------------------------------
G       = ones(n,1)*diag(g)';
i       = 1:size(P.G,2);
G(:,i)  = G(:,i).*exp(P.G);


% Modulatory effects of sp depolarisation on recurrent inhibition
%--------------------------------------------------------------------------
if isfield(P,'M')
    G(:,2) = G(:,2).*exp(-P.M*32*S(:,2));
end

 
% Motion of states: f(x)
%--------------------------------------------------------------------------
 
% Conductance
%==========================================================================
 
% Granular layer (excitatory interneurons): spiny stellate: Hidden causes
%--------------------------------------------------------------------------
u      = G(:,1).*S(:,1) + g(1,3)*S(:,3) + g(1,2)*S(:,2) + A{1}*S(:,2) + U;
f(:,2) = (u - x(:,2)).*T(:,1);
 
% Supra-granular layer (superficial pyramidal cells): Hidden causes - error
%--------------------------------------------------------------------------
u      = g(2,1)*S(:,1) + G(:,2).*S(:,2) + g(2,3)*S(:,3) - A{3}*S(:,4);
f(:,4) = (u - x(:,4)).*T(:,2);
 
% Supra-granular layer (inhibitory interneurons): Hidden states - error
%--------------------------------------------------------------------------
u      = g(3,1)*S(:,1) + g(3,4)*S(:,4) + G(:,3).*S(:,3) + g(3,2)*S(:,2) - A{4}*S(:,4);
f(:,6) = (u - x(:,6)).*T(:,3);
 
% Infra-granular layer (deep pyramidal cells): Hidden states
%--------------------------------------------------------------------------
u      = G(:,4).*S(:,4) + g(4,3)*S(:,3) + A{2}*S(:,2);
f(:,8) = (u - x(:,8)).*T(:,4);
 
% Voltage
%==========================================================================
f(:,1) = x(:,2) - x(:,1).*T(:,1);
f(:,3) = x(:,4) - x(:,3).*T(:,2);
f(:,5) = x(:,6) - x(:,5).*T(:,3);
f(:,7) = x(:,8) - x(:,7).*T(:,4);
f      = spm_vec(f);
 
 
if nargout < 2; return, end
 
 
% delays
%==========================================================================
% Delay differential equations can be integrated efficiently (but
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = Q(d)f(x(t))
%
%    J(d)         = Q(d)df/dx
%--------------------------------------------------------------------------
% Implement: dx(t)/dt = f(x(t - d)) = inv(1 + D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
[Q,J] = spm_dcm_delay(P,M);
 
 
return
 
% notes and alpha function (kernels)
%==========================================================================
% x   = t*exp(k*t)
% x'  = exp(k*t) + k*t*exp(k*t)
%     = exp(k*t) + k*x
% x'' = 2*k*exp(k*t) + k^2*t*exp(k*t)
%     = 2*k*(x' - k*x) + k^2*x
%     = 2*k*x' - k^2*x
