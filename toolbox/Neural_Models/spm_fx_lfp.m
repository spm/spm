function [f,J] = spm_fx_lfp(x,u,P,M)
% state equations for a neural mass model of erps
% FORMAT [f,J] = spm_fx_lfp(x,u,P,M)
% x      - state vector
%   x(:,1)  - voltage (spiny stellate cells)
%   x(:,2)  - voltage (pyramidal cells)         +ve
%   x(:,3)  - voltage (pyramidal cells)         -ve
%   x(:,4)  - current (spiny stellate cells)    +ve 
%   x(:,5)  - current (pyramidal cells)         +ve
%   x(:,6)  - current (pyramidal cells)         -ve
%   x(:,7)  - voltage (inhibitory interneurons) +ve
%   x(:,8)  - current (inhibitory interneurons) +ve
%   x(:,9)  - voltage (pyramidal cells)
%   x(:,10) - voltage (inhibitory interneurons) -ve
%   x(:,11) - current (inhibitory interneurons) -ve
%   x(:,12) - voltage (inhibitory interneurons)
%
%   x(:,13) - slow potassium conductance
%
% f    = dx(t)/dt  = f(x(t))
% J    = df/dx
%
% Fixed parameter scaling [Defaults]
%
%  E = [32 16 4];             % extrinsic rates (forward, backward, lateral)
%  G = [1 1 1/2 1/2 1/8]*128; % intrinsic rates (g1, g2, g3, g4, g5)
%  D = [2 16];                % propagation delays (intrinsic, extrinsic)
%  H = [4 32];                % receptor densities (excitatory, inhibitory)
%  T = [4 16];                % synaptic constants (excitatory, inhibitory)
%  R = [2 1];                 % parameters of static nonlinearity
%
%__________________________________________________________________________
%
% This is a simplified version of spm_fx_erp
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_lfp.m 5369 2013-03-28 20:09:27Z karl $

% check if intrinsic connections are free parameters
%--------------------------------------------------------------------------
try, P.H; catch, P.H = 0; end

% get dimensions and configure state variables
%--------------------------------------------------------------------------
x    = spm_unvec(x,M.x);       % neuronal states
n    = size(x,1);              % number of sources
s    = size(x,2);              % number of states

% [default] fixed parameters
%--------------------------------------------------------------------------
E    = [32 16 4];              % extrinsic rates (forward, backward, lateral)
G    = [1 1 1/2 1/2 1/32]*128; % intrinsic rates (g1, g2 g3, g4)
D    = [2 4];                  % propagation delays (intrinsic, extrinsic)
H    = [8 32];                 % receptor densities (excitatory, inhibitory)
T    = [4 16];                 % synaptic constants (excitatory, inhibitory)
R    = [1 2];                  % parameters of static nonlinearity

% [specified] fixed parameters
%--------------------------------------------------------------------------
if isfield(M,'pF')
    try, E  = M.pF.E; end
    try, G  = M.pF.H; end
    try, D  = M.pF.D; end
    try, H  = M.pF.G; end
    try, T  = M.pF.T; end
    try, R  = M.pF.R; end
end

% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*E(1);
A{2} = exp(P.A{2})*E(2);
A{3} = exp(P.A{3})*E(3);
C    = exp(P.C);
G    = exp(P.H)*diag(G);
 
% intrinsic connectivity and parameters
%--------------------------------------------------------------------------
Te   = T(1)/1000*exp(P.T(:,1));      % excitatory time constants
Ti   = T(2)/1000*exp(P.T(:,2));      % inhibitory time constants
Tk   = 512/1000;                     % slow potassium
He   = H(1)*exp(P.G);                % excitatory receptor density
Hi   = H(2);                         % inhibitory receptor density


% pre-synaptic inputs: s(V) with threshold adaptation
%--------------------------------------------------------------------------
R      = R.*exp(P.R);
x      = x';
X      = x;
X(1,:) = X(1,:) - X(13,:);
S      = 1./(1 + exp(-R(1)*(X - R(2)))) - 1./(1 + exp(R(1)*R(2)));
dSdx   = R(1)*exp(-R(1)*(max(X,-128) - R(2)))./(1 + exp(-R(1)*(X - R(2)))).^2;
 
% input
%==========================================================================
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U = u(:)*32;
    
else
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:);
end


% State: f(x) and Jacobian dfdx
%==========================================================================
 
% NB: activity-dependent reduction in inhibitory effective time-constant
%--------------------------------------------------------------------------
Ti    = 4/1000 + Ti;
 
% intrinsic coupling
%--------------------------------------------------------------------------
for i = 1:n
 
    % synaptic dynamics - dfdx
    %----------------------------------------------------------------------
    ke    = -2/Te(i);
    ki    = -2/Ti(i);
    Ke    = -1/(Te(i)^2);
    Ki    = -1/(Ti(i)^2);
    dfdx  = [0   0   0   1   0   0   0   0   0   0   0   0   0
             0   0   0   0   1   0   0   0   0   0   0   0   0
             0   0   0   0   0   1   0   0   0   0   0   0   0
             Ke  0   0   ke  0   0   0   0   0   0   0   0   0
             0   Ke  0   0   ke  0   0   0   0   0   0   0   0
             0   0   Ki  0   0   ki  0   0   0   0   0   0   0
             0   0   0   0   0   0   0   1   0   0   0   0   0
             0   0   0   0   0   0  Ke   ke  0   0   0   0   0
             0   0   0   0   1   -1  0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   1   0   0
             0   0   0   0   0   0   0   0   0  Ki   ki  0   0
             0   0   0   0   0   0   0   1   0   0  -1   0   0
             0   0   0   0   0   0   0   0   0   0   0   0  -1/Tk];
         
    % intrinsic afferents - dfdS
    %----------------------------------------------------------------------
    Ke    = He(i)/Te(i);
    Ki    = Hi/Ti(i);
    
    j1    = Ke*G(i,1);
    j2    = Ke*G(i,2);
    j3    = Ke*G(i,3);
    j4    = Ki*G(i,4);
    j5    = Ki*G(i,5);
    dfdS  = [0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   j1  0   0   0   0
             j2  0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   j4  0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   j3  0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0  j5   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
            4/Tk 0   0   0   0   0   0   0   0   0   0   0   0];
     
    % motion and Jacobian
    %----------------------------------------------------------------------
    dsdx       = diag(dSdx(:,i));
    dsdx(1,13) = -dsdx(1,1);
    dfdu       = sparse(4,1,Ke,s,1);
    
    F{i}       = dfdx*x(:,i) + dfdS*S(:,i) + dfdu*U(i);
    J{i,i}     = dfdx + dfdS*dsdx;
    
    % extrinsic afferents 
    %----------------------------------------------------------------------
    for j = 1:n, if i ~= j
    
    k1    = Ke*(A{1}(i,j) + A{3}(i,j));
    k2    = Ke*(A{2}(i,j) + A{3}(i,j));
    k3    = Ki*(A{2}(i,j) + A{3}(i,j));
    dfdS  = [0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   k1  0   0   0   0
             0   0   0   0   0   0   0   0   k2  0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   k3  0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0
             0   0   0   0   0   0   0   0   0   0   0   0   0];
         
    % motion and Jacobian
    %----------------------------------------------------------------------
    F{i}   = F{i} + dfdS*S(:,j);
    J{i,j} = dfdS*diag(dSdx(:,j));
    
    end, end
end
 
% construct motion and Jacobian
%--------------------------------------------------------------------------
for i = 1:n
    k      = (1:n:s*n) + (i - 1);
    f(k,1) = F{i};
    for j  = 1:n
        l         = [1:n:s*n] + (j - 1);
        dfdx(k,l) = J{i,j};
    end    
end
 
% extrinsic and intrinsic delays
%--------------------------------------------------------------------------
De = D(2).*exp(P.D)/1000;
Di = D(1).*exp(P.I)/1000;
De = (eye(n,n) - 1).*De;
Di = (eye(s,s) - 1)*Di;
De = kron(ones(s,s),De);
Di = kron(Di,eye(n,n));
 
D  = Di + De;
 
% Implement: dx(t)/dt = f(x(t + d)) = inv(1 - D.*dfdx)*f(x(t))
%--------------------------------------------------------------------------
D  = spm_inv(speye(n*s,n*s) - D.*dfdx);
f  = D*f;
J  = D*dfdx;
 
return
 
% Equations of motion
%==========================================================================
 
% Supragranular layer (inhibitory interneurons): depolarizing current
%--------------------------------------------------------------------------
f(:,7)  = x(:,8);
f(:,8)  = (He.*((A{2} + A{3})*S(:,9) + G(:,3).*S(:,9)) ...
           - 2*x(:,8) - x(:,7)./Te)./Te;
      
% Supragranular layer (inhibitory interneurons): hyperpolarizing current
%--------------------------------------------------------------------------
f(:,10) = x(:,11);
f(:,11) = (Hi*G(:,5).*S(:,12) ...
           - 2*x(:,11) - x(:,10)./Ti)./Ti;
 
% Granular layer (spiny stellate cells): depolarizing current
%--------------------------------------------------------------------------
f(:,1)  = x(:,4);
f(:,4)  = (He.*((A{1} + A{3})*S(:,9) + G(:,1).*S(:,9) + U) ...
           - 2*x(:,4) - x(:,1)./Te)./Te;
       
% Infra-granular layer (pyramidal cells): depolarizing current
%--------------------------------------------------------------------------
f(:,2)  = x(:,5);
f(:,5)  = (He.*((A{2} + A{3})*S(:,9) + G(:,2).*S(:,1)) ...
           - 2*x(:,5) - x(:,2)./Te)./Te;
 
% Infra-granular layer (pyramidal cells): hyperpolarizing current
%--------------------------------------------------------------------------
f(:,3)  = x(:,6);
f(:,6)  = (Hi*G(:,4).*S(:,12) ...
           - 2*x(:,6) - x(:,3)./Ti)./Ti;
 
% Surpa and Infra-granular layer (pyramidal cells): Voltage
%--------------------------------------------------------------------------
f(:,9)  = x(:,5) - x(:,6);
f(:,12) = x(:,8) - x(:,11);
 
% Granular layer (spiny stellate cells): hyperpolarizing current
%--------------------------------------------------------------------------
f(:,13) = (4*S(:,1) - x(:,13))./Tk;
 
% Jacobian for delays (evaluate numerically for simplicity)
%==========================================================================
dfdx = spm_diff('spm_fx_lfp',x,u,P,1,1);
