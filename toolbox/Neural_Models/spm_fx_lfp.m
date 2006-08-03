function [f] = spm_fx_lfp(x,u,P,s)
% state equations for a neural mass model of erps
% FORMAT [f] = spm_fx_lfp(x,u,P)
% FORMAT [f] = spm_fx_lfp(n)
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
% NB: the first state is actually time but this hidden here.
%
% f        - dx(t)/dt  = f(x(t))
%
% Fixed parameter scaling [Defaults]
%
%  E = [32 16 4];             % extrinsic rates (forward, backward, lateral)
%  G = [1 1 1/2 1/2 1/8]*128; % intrinsic rates (g1, g2 g3, g4, g5)
%  D = [2 16];                % propogation delays (intrinsic, extrinsic)
%  H = [4 32];                % receptor densities (excitatory, inhibitory)
%  T = [4 16];                % synaptic constants (excitatory, inhibitory)
%  R = [2 1];                 % parameters of static nonlinearity
%
%__________________________________________________________________________
%
% This is a simplified version of spm_gx_erp
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%

% get dimensions and configure state variables
%--------------------------------------------------------------------------
s     = length(x);                  % number of states
m     = length(u);  				% number of inputs
n     = length(P.A{1});             % number of sources
x     = reshape(x,n,s);             % neuronal states

% effective extrinsic connectivity
%--------------------------------------------------------------------------
for i = 1:m
	  P.A{1} = P.A{1} + u(i)*P.B{i};              % forward  connections
      P.A{2} = P.A{2} + u(i)*P.B{i};              % backward connections
      P.A{3} = P.A{3} + u(i)*P.B{i};              % lateral  connections
end

% [default] fixed parameters
%--------------------------------------------------------------------------

E = [32 16 4];                % extrinsic rates (forward, backward, lateral)
G = [1 1 1/2 1/2 1/8]*128;    % intrinsic rates (g1, g2 g3, g4)
D = [2 16];                   % propogation delays (intrinsic, extrinsic)
H = [4 32];                   % receptor densities (excitatory, inhibitory)
T = [4 16];                   % synaptic constants (excitatory, inhibitory)
R = [2 1];                    % parameters of static nonlinearity

% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*E(1);
A{2} = exp(P.A{2})*E(2);
A{3} = exp(P.A{3})*E(3);
C    = exp(P.C);
G    = exp(P.G)*diag(G);

% intrinsic connectivity and parameters
%--------------------------------------------------------------------------
Te  =  T(1)/1000*exp(P.T(:,1));      % excitatory time constants
Ti   = T(2)/1000*exp(P.T(:,2));      % inhibitory time constants
Tk   = 512/1000;                     % slow potassium
He   = H(1)*exp(P.H);                % excitatory receptor density
Hi   = H(2);                         % inhibitory receptor density

% pre-synaptic inputs: s(V) with threshold adaptation
%--------------------------------------------------------------------------
R      = R.*exp(P.R);
S      = x;
S(:,1) = S(:,1) - x(:,13);
S      = 1./(1 + exp(-R(1)*(S - R(2)))) - 1./(1 + exp(R(1)*R(2)));

% input
%--------------------------------------------------------------------------
U      = C*u;

% State: f(x)
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
       
% Granular layer (spiny stellate cells): hyperpolarizing current
%--------------------------------------------------------------------------
f(:,13) = x(:,1) - x(:,13)./Tk;

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

% vectorise
%--------------------------------------------------------------------------
f  = f(:);

% Delays
%==========================================================================
if nargin > 3, return, end

% Jacobian for delays (evaluate numerically for simplicity)
%--------------------------------------------------------------------------
dfdx = spm_diff('spm_fx_lfp',x,u,P,1,1);

% extrinsic and intrinsic delays
%--------------------------------------------------------------------------
De = D(2).*exp(P.D)/1000;
Di = D(1)/1000;
De = (eye(n,n) - 1)*De;
Di = sparse([4 5 5 6 8 11],[9 1 9 12 9 12],-Di,s,s);
De = kron(ones(s,s),De);
Di = kron(Di,eye(n,n));

D  = Di + De;

% Implement: dx(t)/dt = f(x(t + d)) = inv(1 - D.*dfdx)*f(x(t))
%--------------------------------------------------------------------------
f  = inv(speye(s,s) - D.*dfdx)*f;

return







% Analytic evaluation of Jacobian: J = df(x)/dx
%==========================================================================
dSdx = 1./(1 + exp(-R*x)).^2.*(R*exp(-R*x));
I    = speye(n,n);
J    = kron(sparse(s,s),sparse(n,n));

% changes in voltage with current
%--------------------------------------------------------------------------
S  = sparse([7 2 1 3 9 10 12],[8 5 4 6 5 11 8],1,s,s); J = J + kron(S,I);
S  = sparse([9 12],[6 11],-1,s,s);                     J = J + kron(S,I);

% synaptic kernel
%--------------------------------------------------------------------------
S  = sparse([8 4 5],[8 4 5],2,s,s);     J = J - kron(S,diag(1./Te));
S  = sparse([6 11],[6 11],2,s,s);       J = J - kron(S,diag(1./Ti));
S  = sparse([8 4 5],[7 1 2],1,s,s);     J = J - kron(S,diag(1./(Te.*Te)));
S  = sparse([6 11],[3 10],1,s,s);       J = J - kron(S,diag(1./(Ti.*Ti)));

% Supragranular layer (inhibitory interneurons)
%--------------------------------------------------------------------------
E  = (A{2} + A{3})*diag(dSdx(:,9)) + diag(G(:,3).*dSdx(:,9));
E  = diag(He./Te)*E;
S  = sparse(8,9,1,s,s);                 J = J + kron(S,E);

% Granular layer (spiny stellate cells)
%--------------------------------------------------------------------------
E  = (A{1} + A{3})*diag(dSdx(:,9)) + diag(G(:,1).*dSdx(:,9));
E  = diag(He./Te)*E;
S  = sparse(4,9,1,s,s);                 J = J + kron(S,E);

% Infraranular layer (pyramidal cells)
%--------------------------------------------------------------------------
E  = (A{2} + A{3})*diag(dSdx(:,9));
E  = diag(He./Te)*E;
S  = sparse(5,9,1,s,s);                 J = J + kron(S,E);

E  = diag(G(:,2).*dSdx(:,1));
E  = diag(He./Te)*E;
S  = sparse(5,1,1,s,s);                 J = J + kron(S,E);

% Infra-granular layer (pyramidal cells)
%--------------------------------------------------------------------------
E  = diag(G(:,4).*dSdx(:,12));
E  = diag(Hi./Ti)*E;
S  = sparse(6,12,1,s,s);                J = J + kron(S,E);

% Suprs-granular layer (inhibitory cells)
%--------------------------------------------------------------------------
E  = diag(G(:,5).*dSdx(:,12));
E  = diag(Hi./Ti)*E;
S  = sparse(11,12,1,s,s);               J = J + kron(S,E);

dfdx = J;


