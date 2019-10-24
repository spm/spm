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
% $Id: spm_fx_cmc_tfm.m 7679 2019-10-24 15:54:07Z spm $
 
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
x  = spm_unvec(x,M.x);                % neuronal states
n  = size(x,1);                       % number of sources


% [default] fixed parameters
%--------------------------------------------------------------------------
E  = [2 1 1 1]*512;                   % extrinsic (forward and backward)  
T  = [16 8 1 2]*16;                   % synaptic rate constants
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
ig   = 2*(1:4);                       % indices of conductance
iv   = ig - 1;                        % indices of voltage
V    = x(:,iv);                       % Voltage
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
% forward  A{1}  sp -> ss
% forward  A{2}  sp -> dp
% backward A{3}  dp -> sp
% backward A{4}  dp -> ii
%--------------------------------------------------------------------------

% Neuronal states (deviations from baseline)
%--------------------------------------------------------------------------
%   x(:,1) - voltage     (spiny stellate cells)
%   x(:,2) - conductance (spiny stellate cells)but
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - current     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%   x(:,7) - voltage     (deep pyramidal cells)
%   x(:,8) - conductance (deep pyramidal cells)
%--------------------------------------------------------------------------
%     ss sp ii dp   % intrinsic connections
%--------------------------------------------------------------------------
g  = [-8 -4  -4  0;  % ss
       4 -8  -2  0;  % sp
       4  2  -4 -2;  % ii
       0  1  -2 -4]; % dp

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

% Afferents
%==========================================================================
u      = zeros(n,4);             % intrinsic - inhibitory
v      = zeros(n,4);             % intrinsic - excitatory
w      = zeros(n,4);             % extrinsic - excitatory
 
% Granular layer (excitatory interneurons): spiny stellate: Hidden causes
%--------------------------------------------------------------------------
u(:,1) = G(:,1).*S(:,1) + g(1,3)*S(:,3) + g(1,2)*S(:,2);
w(:,1) = A{1}*S(:,2) + U;
 
% Supra-granular layer (superficial pyramidal cells): Hidden causes - error
%--------------------------------------------------------------------------
u(:,2) = G(:,2).*S(:,2) + g(2,3)*S(:,3);
v(:,2) = g(2,1)*S(:,1);
w(:,2) = A{3}*S(:,4);

% Supra-granular layer (inhibitory interneurons): Hidden states - error
%--------------------------------------------------------------------------
u(:,3) = G(:,3).*S(:,3);
v(:,3) = g(3,1)*S(:,1) + g(3,4)*S(:,4) + g(3,2)*S(:,2);
w(:,3) = A{4}*S(:,4);

% Infra-granular layer (deep pyramidal cells): Hidden states
%--------------------------------------------------------------------------
u(:,4) = G(:,4).*S(:,4) + g(4,3)*S(:,3);
v(:,4) = g(4,2)*S(:,2);
w(:,4) = A{2}*S(:,2);


if nargin > 4; f = u; J = v; Q = w; return, end
 
% Conductance and voltage
%==========================================================================
f(:,ig) = (u + v + w - x(:,ig)).*T;
f(:,iv) = x(:,ig) - x(:,iv).*T;
f       = spm_vec(f);
 

if nargout < 2; return, end

% evaluate Jacobian
%==========================================================================

% derivatives of firing with respect to voltage
%--------------------------------------------------------------------------
dS     = (R*exp(B - R*V))./(exp(B - R*V) + 1).^2;

% intrinsic connectivity
%--------------------------------------------------------------------------
for i = 1:4
    for j = 1:4
        if i == j
            J{ig(i),iv(i)} = diag(G(:,i).*dS(:,j).*T(:,i));
        else
            J{ig(i),iv(j)} = diag(g(i,j) *dS(:,j).*T(:,i));
        end
    end
end

% extrinsic connectivity
%--------------------------------------------------------------------------
J{2,3} = A{1}*diag(dS(:,2).*T(:,1)) + J{2,3};
J{4,7} = A{3}*diag(dS(:,4).*T(:,2)) + J{4,7};
J{6,7} = A{4}*diag(dS(:,4).*T(:,3)) + J{6,7};
J{8,3} = A{2}*diag(dS(:,2).*T(:,4)) + J{8,3};

% conductance decay
%--------------------------------------------------------------------------
for i = 1:4
    J{ig(i),ig(i)} = -diag(T(:,i));
    J{iv(i),iv(i)} = -diag(T(:,i));
    J{iv(i),ig(i)} =  speye(n,n);
end
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






