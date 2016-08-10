function [f,J,Q] = spm_fx_cmc_tfm(x,u,P,M,OPT)
% state equations for a neural mass model (canonical microcircuit)
% FORMAT [f,J,D]  = spm_fx_cmc_tfm(x,u,P,M)
% FORMAT [f,J]    = spm_fx_cmc_tfm(x,u,P,M)
% FORMAT [f]      = spm_fx_cmc_tfm(x,u,P,M)
%
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
% f  - dx(t)/dt  = f(x(t))
% J  - df(t)/dx(t)
% D  - delay operator dx(t)/dt = f(x(t - d))
%                                    = D(d)*f(x(t))
%
% FORMAT [u,v,w] = spm_fx_cmc_tfm(x,u,P,M,'activity')
% u  - intrinsic presynaptic input (inhibitory)
% v  - intrinsic presynaptic input (excitatory)
% w  - extrinsic presynaptic input
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
% $Id: spm_fx_cmc_tfm.m 6856 2016-08-10 17:55:05Z karl $
 
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
x  = spm_unvec(x,M.x);                % neuronal states
n  = size(x,1);                       % number of sources


% [default] fixed parameters
%--------------------------------------------------------------------------
E  = [2 1 1 1]*512;                   % extrinsic (forward and backward)  
T  = [256 128 16 32];                    % synaptic rate constants
R  = 1;                               % gain of activation function
B  = 0;                               % baseline firing

% Extrinsic connections
%--------------------------------------------------------------------------
% ss = spiny stellate
% sp = superficial pyramidal
% dp = deep pyramidal
% ii = inhibitory interneurons
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*E(1);              % forward  connections (sp -> ss)
A{2} = exp(P.A{2})*E(2);              % forward  connections (sp -> dp)
A{3} = exp(P.A{3})*E(3);              % backward connections (dp -> sp)
A{4} = exp(P.A{4})*E(4);              % backward connections (dp -> ii)

% input connections
%--------------------------------------------------------------------------
C    = exp(P.C);
 
% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
V    = x(:,1:2:end);                  % Voltage
F    = 1./(1 + exp(B - R*V));         % firing rate
S    = F - 1/(1 + exp(B));            % deviation 

% input
%==========================================================================
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U   = u(:)*128;
    M.m = size(U,1);
    
else
    % exogenous input
    %----------------------------------------------------------------------
    U   = C*u(:)*2;
    M.m = size(C,2);
    
end
clear u

 
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
g  = [-8 -4  -4  0;  % ss
       4 -8  -2  0;  % sp
       4  2  -4 -2;  % ii
       0  2  -2 -4]; % dp

g  = g*256*exp(P.S);

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
u     = zeros(n,4);             % intrinsic – inhibitory
v     = zeros(n,4);             % intrinsic – excitatory
w     = zeros(n,4);             % extrinsic – excitatory

% Conductance
%==========================================================================
 
% Granular layer (excitatory interneurons): spiny stellate: Hidden causes
%--------------------------------------------------------------------------
u(:,1) = G(:,1).*S(:,1) + g(1,3)*S(:,3) + g(1,2)*S(:,2);
w(:,1) = A{1}*S(:,2) + U;
f(:,2) = (u(:,1) + w(:,1) - x(:,2)).*T(:,1);
 
% Supra-granular layer (superficial pyramidal cells): Hidden causes - error
%--------------------------------------------------------------------------
u(:,2) = G(:,2).*S(:,2) + g(2,3)*S(:,3);
v(:,2) = g(2,1)*S(:,1);
w(:,2) = A{3}*S(:,4);
f(:,4) = (u(:,2) + v(:,2) + w(:,2) - x(:,4)).*T(:,2);
 
% Supra-granular layer (inhibitory interneurons): Hidden states - error
%--------------------------------------------------------------------------
u(:,3) = G(:,3).*S(:,3);
v(:,3) = g(3,1)*S(:,1) + g(3,4)*S(:,4) + g(3,2)*S(:,2);
w(:,3) = A{4}*S(:,4);
f(:,6) = (u(:,3) + v(:,3) + w(:,3) - x(:,6)).*T(:,3);
 
% Infra-granular layer (deep pyramidal cells): Hidden states
%--------------------------------------------------------------------------
u(:,4) = G(:,4).*S(:,4) + g(4,3)*S(:,3);
v(:,4) = g(4,2)*S(:,2);
w(:,4) = A{2}*S(:,2);
f(:,8) = (u(:,4) + w(:,4) - x(:,8)).*T(:,4);

if nargin > 4; f = -u; J = v; Q = w; return, end
 
% Voltage
%==========================================================================
f(:,1) = x(:,2) - x(:,1).*T(:,1);
f(:,3) = x(:,4) - x(:,3).*T(:,2);
f(:,5) = x(:,6) - x(:,5).*T(:,3);
f(:,7) = x(:,8) - x(:,7).*T(:,4);
f      = spm_vec(f);
 

if nargout < 2; return, end

% evaluate Jacobian
%==========================================================================

% derivatives of firing with respect to voltage
%--------------------------------------------------------------------------
dS     = (R*exp(B - R*V))./(exp(B - R*V) + 1).^2;

% intrinsic connectivity
%--------------------------------------------------------------------------
J{2,3} = diag(g(1,2) *dS(:,2).*T(:,1));
J{2,5} = diag(g(1,3) *dS(:,3).*T(:,1));
J{2,1} = diag(G(:,1).*dS(:,1).*T(:,1));

J{4,5} = diag(g(2,3) *dS(:,3).*T(:,2));
J{4,1} = diag(g(2,1) *dS(:,1).*T(:,2));
J{4,3} = diag(G(:,2).*dS(:,2).*T(:,2));

J{6,3} = diag(g(3,2) *dS(:,2).*T(:,3));
J{6,7} = diag(g(3,4) *dS(:,4).*T(:,3));
J{6,1} = diag(g(3,1) *dS(:,1).*T(:,3));
J{6,5} = diag(G(:,3).*dS(:,3).*T(:,3));

J{8,3} = diag(g(4,2) *dS(:,2).*T(:,4));
J{8,5} = diag(g(4,3) *dS(:,3).*T(:,4));
J{8,7} = diag(G(:,4).*dS(:,4).*T(:,4));

% extrinsic connectivity
%--------------------------------------------------------------------------
J{2,3} = A{1}*diag(dS(:,2).*T(:,1)) + J{2,3};
J{4,7} = A{3}*diag(dS(:,4).*T(:,2));
J{6,7} = A{4}*diag(dS(:,4).*T(:,3)) + J{6,7};
J{8,3} = A{2}*diag(dS(:,2).*T(:,4));

% conductance decay
%--------------------------------------------------------------------------
J{2,2} = -diag(T(:,1));
J{4,4} = -diag(T(:,2));
J{6,6} = -diag(T(:,3));
J{8,8} = -diag(T(:,4));

% voltage decay
%--------------------------------------------------------------------------
J{1,1} = -diag(T(:,1)); J{1,2} = speye(n,n);
J{3,3} = -diag(T(:,2)); J{3,4} = speye(n,n);
J{5,5} = -diag(T(:,3)); J{5,6} = speye(n,n);
J{7,7} = -diag(T(:,4)); J{7,8} = speye(n,n);

J     = spm_cat(J);

if nargout < 3; return, end
 
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
[Q,J] = spm_dcm_delay(P,M,J,0);
 
 
return
 
% notes and alpha function (kernels)
%==========================================================================
% x   = t*exp(k*t)
% x'  = exp(k*t) + k*t*exp(k*t)
%     = exp(k*t) + k*x
% x'' = 2*k*exp(k*t) + k^2*t*exp(k*t)
%     = 2*k*(x' - k*x) + k^2*x
%     = 2*k*x' - k^2*x



% notes on Eigensolutions
%==========================================================================
clear J
syms g11 g22 g33 g44 g12 g23 T1 T2 T3 T4

% rate constants
%--------------------------------------------------------------------------
T  = [T1, T2, T3, T4];

% intrinsic connectivity
%--------------------------------------------------------------------------
g  = [-g11  -g12  -g12  0;    % ss
       g12  -g22  -g23  0;    % sp
       g12   g23  -g33  g23;  % ii
       0     0    -g23 -g44]; % dp

% intrinsic connectivity
%--------------------------------------------------------------------------
J = [[-diag(T)    eye(4) ]
     [ diag(T)*g -diag(T)]];

% i = [1 3 5 7 2 4 6 8]; J(i,i) = 
% 
% [     -T1,       0,       0,       0,   1,   0,   0,   0]
% [       0,     -T2,       0,       0,   0,   1,   0,   0]
% [       0,       0,     -T3,       0,   0,   0,   1,   0]
% [       0,       0,       0,     -T4,   0,   0,   0,   1]
% [ -T1*g11, -T1*g12, -T1*g12,       0, -T1,   0,   0,   0]
% [  T2*g12, -T2*g22, -T2*g23,       0,   0, -T2,   0,   0]
% [  T3*g12,  T3*g23, -T3*g33,  T3*g23,   0,   0, -T3,   0]
% [       0,       0, -T4*g23, -T4*g44,   0,   0,   0, -T4]






