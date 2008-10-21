function [f,J,Q] = spm_fx_erp(x,u,P,M)
% state equations for a neural mass model of erps
% FORMAT [f,J,D] = spm_fx_erp(x,u,P,M)
% FORMAT [f,J]   = spm_fx_erp(x,u,P,M)
% FORMAT [f]     = spm_fx_erp(x,u,P,M)
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
% J        - df(t)/dx(t)
% D        - delay operator dx(t)/dt = f(x(t - d))
%                                    = D(d)*f(x(t))
%
% Prior fixed parameter scaling [Defaults]
%
%  M.pF.E = [32 16 4];           % extrinsic rates (forward, backward, lateral)
%  M.pF.G = [1 4/5 1/4 1/4]*128; % intrinsic rates (g1, g2 g3, g4)
%  M.pF.D = [2 16];              % propogation delays (intrinsic, extrinsic)
%  M.pF.H = [4 32];              % receptor densities (excitatory, inhibitory)
%  M.pF.T = [8 16];              % synaptic constants (excitatory, inhibitory)
%  M.pF.R = 0.56;                % parameter of static nonlinearity
%
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_erp.m 2374 2008-10-21 18:52:29Z karl $




% get dimensions and configure state variables
%--------------------------------------------------------------------------
n  = length(P.A{1});             % number of sources
x  = spm_unvec(x,M.x);           % neuronal states

% [default] fixed parameters
%--------------------------------------------------------------------------
try
    E  = M.pF.E;                % extrinsic rates (forward, backward, lateral)
    G  = M.pF.G;                % intrinsic rates (g1, g2 g3, g4)
    D  = M.pF.D;                % propogation delays (intrinsic, extrinsic)
    H  = M.pF.H;                % receptor densities (excitatory, inhibitory)
    T  = M.pF.T;                % synaptic constants (excitatory, inhibitory)
    R  = M.pF.R;                % parameters of static nonlinearity
catch
    E = [32 16 4];              % extrinsic rates (forward, backward, lateral)
    G = [1 4/5 1/4 1/4]*128;    % intrinsic rates (g1, g2 g3, g4)
    D = [2 32];                 % propogation delays (intrinsic, extrinsic)
    H = [4 32];                 % receptor densities (excitatory, inhibitory)
    T = [8 16];                 % synaptic constants (excitatory, inhibitory)
    R = [2 1]/3;                % parameters of static nonlinearity
end


% test for free parameters on intrinsic connections
%--------------------------------------------------------------------------
G     = ones(n,1)*G;
try
    G = G.*exp(P.G);
end

% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1}  = exp(P.A{1})*E(1);
A{2}  = exp(P.A{2})*E(2);
A{3}  = exp(P.A{3})*E(3);
C     = exp(P.C);

% intrinsic connectivity and parameters
%--------------------------------------------------------------------------
Te    = T(1)/1000*exp(P.T);           % excitatory time constants
Ti    = T(2)/1000;                    % inhibitory time constants
Hi    = H(2);                         % inhibitory receptor density
He    = H(1)*exp(P.H);                % excitatory receptor density

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R     = R.*exp(P.S);
S     = 1./(1 + exp(-R(1)*(x - R(2)))) - 1./(1 + exp(R(1)*R(2)));
dSdx  = 1./(1 + exp(-R(1)*(x - R(2)))).^2.*(R(1)*exp(-R(1)*(x - R(2))));

% exogenous input
%==========================================================================
U     = C*u(:);


% State: f(x)
%==========================================================================

% Supragranular layer (inhibitory interneurons): Voltage & depolarizing current
%--------------------------------------------------------------------------
f(:,7) = x(:,8);
f(:,8) = (He.*((A{2} + A{3})*S(:,9) + G(:,3).*S(:,9)) - 2*x(:,8) - x(:,7)./Te)./Te;

% Granular layer (spiny stellate cells): Voltage & depolarizing current
%--------------------------------------------------------------------------
f(:,1) = x(:,4);
f(:,4) = (He.*((A{1} + A{3})*S(:,9) + G(:,1).*S(:,9) + U) - 2*x(:,4) - x(:,1)./Te)./Te;

% Infra-granular layer (pyramidal cells): depolarizing current
%--------------------------------------------------------------------------
f(:,2) = x(:,5);
f(:,5) = (He.*((A{2} + A{3})*S(:,9) + G(:,2).*S(:,1)) - 2*x(:,5) - x(:,2)./Te)./Te;

% Infra-granular layer (pyramidal cells): hyperpolarizing current
%--------------------------------------------------------------------------
f(:,3) = x(:,6);
f(:,6) = (Hi*G(:,4).*S(:,7) - 2*x(:,6) - x(:,3)/Ti)/Ti;

% Infra-granular layer (pyramidal cells): Voltage
%--------------------------------------------------------------------------
f(:,9) = x(:,5) - x(:,6);

f      = f(:);

if nargout == 1; return, end



% Jacobian: J = df(x)/dx
%===========================================================================
I  = speye(n,n);
J  = kron(sparse(9,9),sparse(n,n));

% changes in voltage with current
%--------------------------------------------------------------------------
S  = sparse([7 2 1 3 9],[8 5 4 6 5],1,9,9);  J = J + kron(S,I);
S  = sparse(9,6,-1,9,9);                     J = J + kron(S,I);

% synaptic kernel
%--------------------------------------------------------------------------
S  = sparse([8 4 5],[8 4 5],1,9,9); J = J - kron(S,diag(2./Te));
S  = sparse(6,6,1,9,9);             J = J - kron(S,I)*2/Ti;
S  = sparse([8 4 5],[7 1 2],1,9,9); J = J - kron(S,diag(1./(Te.*Te)));
S  = sparse(6,3,1,9,9);             J = J - kron(S,I)/(Ti*Ti);

% Supragranular layer (inhibitory interneurons)
%--------------------------------------------------------------------------
E  = (A{2} + A{3})*diag(dSdx(:,9)) + diag(G(:,3).*dSdx(:,9));
E  = diag(He./Te)*E;
S  = sparse(8,9,1,9,9); J = J + kron(S,E);

% Granular layer (spiny stellate cells)
%--------------------------------------------------------------------------
E  = (A{1} + A{3})*diag(dSdx(:,9)) + diag(G(:,1).*dSdx(:,9));
E  = diag(He./Te)*E;
S  = sparse(4,9,1,9,9); J = J + kron(S,E);

% Infra-granular layer (pyramidal cells)
%--------------------------------------------------------------------------
E  = (A{2} + A{3})*diag(dSdx(:,9));
E  = diag(He./Te)*E;
S  = sparse(5,9,1,9,9); J = J + kron(S,E);

E  = diag(G(:,2).*dSdx(:,1));
E  = diag(He./Te)*E;
S  = sparse(5,1,1,9,9); J = J + kron(S,E);

% Infra-granular layer (pyramidal cells)
%--------------------------------------------------------------------------
E  = diag(G(:,4).*dSdx(:,7));
E  = Hi*E/Ti;
S  = sparse(6,7,1,9,9); J = J + kron(S,E);

if nargout == 2; return, end




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


De = D(2).*exp(P.D)/1000;
Di = D(1)/1000;
De = (1 - speye(n,n)).*De;
Di = (1 - speye(9,9)).*Di;
De = kron(ones(9,9),De);
Di = kron(Di,speye(n,n));
D  = Di + De;

% Implement: dx(t)/dt = f(x(t - d)) = inv(1 + D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = inv(speye(length(J)) + D.*J);





