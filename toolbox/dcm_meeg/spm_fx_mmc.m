function [f,J,Q] = spm_fx_mmc(x,u,P,M)
% state equations for a neural mass model of motor cortex
% Bhatt et al. 2016 Neuroimage
%
% FORMAT [f,J,D] = spm_fx_mmc(x,u,P,M)
% FORMAT [f,J]   = spm_fx_mmc(x,u,P,M)
% FORMAT [f]     = spm_fx_mmc(x,u,P,M)
% x      - state vector
%   x(:,1) - voltage     (middle pyramidal cells)
%   x(:,2) - conductance (middle pyramdidal cells)
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
% T  = synaptic time constants
% S  = slope of sigmoid activation function
%
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging
 
% Bernadette van Wijk
% $Id: spm_fx_mmc.m 7185 2017-10-11 10:10:04Z spm $


% get dimensions and configure state variables
%--------------------------------------------------------------------------
x  = spm_unvec(x,M.x);            % neuronal states
n  = size(x,1);                   % number of sources


% [default] fixed parameters
%--------------------------------------------------------------------------
G  = [2 4 2 2 2 2 2 2 2 2 4 2 2 2]*200;         % intrinsic connections
T  = [3 2 12 18];                               % synaptic time constants [mp sp ii dp]

if isfield(M,'MMC_G'); G = M.MMC_G; end
if isfield(M,'MMC_T'); T = M.MMC_T; end

 
% [specified] fixed parameters
%--------------------------------------------------------------------------
if isfield(M,'pF')
    try, E = M.pF.E; end
    try, G = M.pF.G; end
    try, T = M.pF.T; end
    try, R = M.pF.R; end
end
 
 
% Extrinsic connections
%--------------------------------------------------------------------------
% mp = middle pyramidal
% sp = superficial pyramidal
% dp = deep pyramidal
% ii = inhibitory interneurons
%--------------------------------------------------------------------------
if n > 1
    E    = [1 1/2 1 1/2]*200;     % extrinsic (forward and backward)
    A{1} = exp(P.A{1})*E(1);      % forward  connections (sp -> mp)
    A{2} = exp(P.A{2})*E(2);      % forward  connections (sp -> sp)
    A{3} = exp(P.A{3})*E(3);      % backward connections (dp -> dp)
    
else
    A    = {0,0,0,0};
end

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
R    = 2/3;                      % slope of sigmoid activation function
B    = 0;                        % bias or background (sigmoid)
R    = R.*exp(P.S);              % gain of activation function
F    = 1./(1 + exp(-R*x + B));   % firing rate
S    = F - 1/(1 + exp(B));       % deviation from baseline firing

% input
%==========================================================================
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U = u(:)*1536;
    
else
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:)*32;
    
end

 
% time constants and intrinsic connections
%==========================================================================
T    = ones(n,1)*T/1000;
G    = ones(n,1)*G;

% extrinsic connections
%--------------------------------------------------------------------------
% forward  (i)   2  sp -> mp (+ve)
% forward  (ii)  1  sp -> dp (+ve)
% backward (i)   2  dp -> sp (-ve)
% backward (ii)  1  dp -> ii (-ve)
%--------------------------------------------------------------------------
% free parameters on time constants and intrinsic connections
%--------------------------------------------------------------------------
% G(:,1)  mp -> mp (-ve self)  4
% G(:,2)  mp -> sp (+ve rec )  4
% G(:,3)  ii -> mp (-ve rec )  4
% G(:,4)  ii -> ii (-ve self)  4
% G(:,5)  mp -> ii (+ve rec )  4
% G(:,6)  dp -> ii (+ve rec )  2
% G(:,7)  sp -> sp (-ve self)  4
% G(:,8)  sp -> mp (+ve rec )  4
% G(:,9)  ii -> dp (-ve rec )  2
% G(:,10) dp -> dp (-ve self)  1
% G(:,11) sp -> dp (+ve rec)  2
% G(:,12) ii -> sp (-ve rec)  4
% G(:,13) sp -> ii (+ve rec)  4
% G(:,14) dp -> sp (+ve rec)  2
%--------------------------------------------------------------------------
% Neuronal states (deviations from baseline firing)
%--------------------------------------------------------------------------
%   S(:,1) - voltage     (middle pyramidal cells)
%   S(:,2) - conductance (middle pyramidal cells)
%   S(:,3) - voltage     (superficial pyramidal cells)
%   S(:,4) - conductance (superficial pyramidal cells)
%   S(:,5) - current     (inhibitory interneurons)
%   S(:,6) - conductance (inhibitory interneurons)
%   S(:,7) - voltage     (deep pyramidal cells)
%   S(:,8) - conductance (deep pyramidal cells)
%--------------------------------------------------------------------------
j     = [1 2 3 4];
for i = 1:size(P.T,2)
    T(:,j(i)) = T(:,j(i)).*exp(P.T(:,i));
end

j = 1:14;
for i = 1:size(P.G,2)
    G(:,j(i)) = G(:,j(i)).*exp(P.G(:,i));
end

% Modulatory effects of dp depolarisation on intrinsic connection j(1)
%--------------------------------------------------------------------------
if isfield(P,'M')
    G(:,j(1)) = G(:,j(1)).*exp(-P.M*32*S(:,7));
end

 
% Motion of states: f(x)
%--------------------------------------------------------------------------
 
% Conductance
%==========================================================================
 
% Middle layer (middle pyramidal cells): Hidden causes
%--------------------------------------------------------------------------
u      =   A{1}*S(:,3)+ U;
u      = - G(:,1).*S(:,1) - G(:,3).*S(:,5) + G(:,8).*S(:,3) + u;
f(:,2) =  (u - 2*x(:,2) - x(:,1)./T(:,1))./T(:,1);
 
% Supra-granular layer (superficial pyramidal cells): Hidden causes - error
%--------------------------------------------------------------------------
u      =  A{2}*S(:,3);% + U;
u      =  - G(:,7).*S(:,3) + G(:,2).*S(:,1) - G(:,12).*S(:,5) + G(:,14).*S(:,7) + u;
f(:,4) =  (u - 2*x(:,4) - x(:,3)./T(:,2))./T(:,2);
 
% Supra-granular layer (inhibitory interneurons): Hidden states - error
%--------------------------------------------------------------------------
u      =  0; %- A{4}*S(:,7);
u      =  - G(:,4).*S(:,5) + G(:,5).*S(:,1) + G(:,6).*S(:,7) + G(:,13).*S(:,3) + u;
f(:,6) =  (u - 2*x(:,6) - x(:,5)./T(:,3))./T(:,3);
 
% Infra-granular layer (deep pyramidal cells): Hidden states
%--------------------------------------------------------------------------
u      =   A{3}*S(:,7);% + U;
u      = - G(:,10).*S(:,7) - G(:,9).*S(:,5) + G(:,11).*S(:,3) + u;
f(:,8) =  (u - 2*x(:,8) - x(:,7)./T(:,4))./T(:,4);
 
% Voltage
%==========================================================================
f(:,1) = x(:,2);
f(:,3) = x(:,4);
f(:,5) = x(:,6);
f(:,7) = x(:,8);
f      = spm_vec(f);
 


if nargout < 2; return, end

% Jacobian
%==========================================================================
if isfield(M,'x'), x = spm_vec(M.x); else,  x = sparse(M.n,1); end
if isfield(M,'u'), u = spm_vec(M.u); else,  u = sparse(M.m,1); end
J  = spm_diff(M.f,x,u,P,M,1);


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
Q  = spm_dcm_delay(P,M,J);
 

return
 
% notes and alpha function (kernels)
%==========================================================================
% x   = t*exp(k*t)
% x'  = exp(k*t) + k*t*exp(k*t)
%     = exp(k*t) + k*x
% x'' = 2*k*exp(k*t) + k^2*t*exp(k*t)
%     = 2*k*(x' - k*x) + k^2*x
%     = 2*k*x' - k^2*x
