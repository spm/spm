function [f,J,Q] = spm_fx_cmc(x,u,P,M)
% state equations for a neural mass model (canonical microcircuit)
% FORMAT [f,J,D] = spm_fx_cmc(x,u,P,M)
% FORMAT [f,J]   = spm_fx_cmc(x,u,P,M)
% FORMAT [f]     = spm_fx_cmc(x,u,P,M)
% x      - state vector
%   x(:,1) - voltage     (spiny stellate cells)
%   x(:,2) - conductance (spiny stellate cells)
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - current     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%   x(:,7) - voltage     (deep pyramidal cells)
%   x(:,8) - conductance (deep pyramidal cells)
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
% D        - delay operator dx(t)/dt = f(x(t - d))
%                                    = D(d)*f(x(t))
%
% Prior fixed parameter scaling [Defaults]
%
% E  = (forward, backward, lateral) extrinsic rates 
% G  = intrinsic rates
% D  = propagation delays (intrinsic, extrinsic)
% H  = overall synaptic kinetics
% T  = synaptic time constants
% R  = slope of sigmoid activation function
%
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_cmc.m 4348 2011-06-10 20:50:23Z karl $
 
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
M.u   = u;                          % place inputs in M
x     = spm_unvec(x,M.x);           % neuronal states
[n m] = size(x);                    % number of sources and states  
 
 
% [default] fixed parameters
%--------------------------------------------------------------------------
E  = [1 1/2 1 1/2]*200;             % extrinsic (forward and backward)  
G  = [4 4 4 4 4 2 4 4 2 1]*200;     % intrinsic connections
D  = [1 16];                        % delays (intrinsic, extrinsic)
T  = [2 2 16 28];                   % synaptic time constants
R  = 2/3;                           % slope of sigmoid activation function
 

% [specified] fixed parameters
%--------------------------------------------------------------------------
try, E = M.pF.E; end
try, G = M.pF.G; end
try, D = M.pF.D; end
try, T = M.pF.T; end
try, R = M.pF.R; end
 
 
% Extrinsic connections
%--------------------------------------------------------------------------
% ss = spiny stellate
% sp = superficial pyramidal
% dp = deep pyramidal
% ii = inhibitory interneurons
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*E(1);           % forward  connections (sp -> ss)
A{2} = exp(P.A{2})*E(2);           % forward  connections (sp -> dp)
A{3} = exp(P.A{3})*E(3);           % backward connections (dp -> sp)
A{4} = exp(P.A{4})*E(4);           % backward connections (dp -> ii)
C    = exp(P.C);
 
% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R    = R.*exp(P.S);
S    = 1./(1 + exp(-R*x)) - 1/2;
 
% exogenous input
%--------------------------------------------------------------------------
U    = C*u(:);
 
 
% time constants and intrinsic connections
%==========================================================================
T    = ones(n,1)*T/1000;
G    = ones(n,1)*G;
 
% free parameters on time constants and intrinsic connections
%--------------------------------------------------------------------------
% G(:,1)  ss -> ss
% G(:,2)  sp -> ss
% G(:,3)  ii -> ss
% G(:,4)  ii -> ii
% G(:,5)  ss -> ii
% G(:,6)  dp -> ii
% G(:,7)  sp -> sp
% G(:,8)  ss -> sp
% G(:,9)  ii -> dp
% G(:,10) dp -> dp
%--------------------------------------------------------------------------
for i = 1:size(P.T,2)
    T(:,i) = T(:,i).*exp(P.T(:,i));
end
for i = 1:size(P.G,2)
    G(:,i) = G(:,i).*exp(P.G(:,i));
end
 
 
% Motion of states: f(x)
%--------------------------------------------------------------------------
 
% Conductance
%==========================================================================
 
% Granular layer (excitatory interneurons): spiny stellate: Hidden causes
%--------------------------------------------------------------------------
u      =   A{1}*S(:,3) + U;
u      = - G(:,1).*S(:,1) - G(:,3).*S(:,5) - G(:,2).*S(:,3) + u;
f(:,2) = (u - 2*x(:,2) - x(:,1)./T(:,1))./T(:,1);
 
% Supra-granular layer (superficial pyramidal cells): Hidden causes - error
%--------------------------------------------------------------------------
u      = - A{3}*S(:,7);
u      =   G(:,8).*S(:,1) - G(:,7).*S(:,3) + u;
f(:,4) = (u - 2*x(:,4) - x(:,3)./T(:,2))./T(:,2);
 
% Supra-granular layer (inhibitory interneurons): Hidden states - error
%--------------------------------------------------------------------------
u      = - A{4}*S(:,7);
u      =   G(:,5).*S(:,1) + G(:,6).*S(:,7) - G(:,4).*S(:,5) + u;
f(:,6) = (u - 2*x(:,6) - x(:,5)./T(:,3))./T(:,3);
 
% Infra-granular layer (deep pyramidal cells): Hidden states
%--------------------------------------------------------------------------
u      =   A{2}*S(:,3);
u      = - G(:,10).*S(:,7) - G(:,9).*S(:,5) + u;
f(:,8) = (u - 2*x(:,8) - x(:,7)./T(:,4))./T(:,4);
 
% Voltage
%==========================================================================
f(:,1) = x(:,2);
f(:,3) = x(:,4);
f(:,5) = x(:,6);
f(:,7) = x(:,8);
f      = spm_vec(f);
 
 
if nargout == 1; return, end
 
 
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
De = exp(P.D);
Di = diag(diag(De));
De = De - Di;
De = De*D(2)/1000;
Di = Di*D(1)/1000;
De = kron(ones(m,m),De);
Di = kron(ones(m,m) - speye(m,m),Di);
D  = Di + De;
 
% Implement: dx(t)/dt = f(x(t - d)) = inv(1 + D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
[Q,J] = spm_dcm_delay(M,P,D);
 
 
return
 
% notes and alpha function (kernels)
%==========================================================================
% x   = t*exp(k*t)
% x'  = exp(k*t) + k*t*exp(k*t)
%     = exp(k*t) + k*x
% x'' = 2*k*exp(k*t) + k^2*t*exp(k*t)
%     = 2*k*(x' - k*x) + k^2*x
%     = 2*k*x' - k^2*x
