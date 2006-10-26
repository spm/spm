function [f,J,D] = spm_fx_erp(x,u,P)
% state equations for a neural mass model of erps
% FORMAT [f,J,D] = spm_fx_erp(x,u,P)
% FORMAT [f,D*J] = spm_fx_erp(x,u,P)
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
% NB: the first state is actually time but this hidden here.
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
% D        - delay operator dx(t)/dt = f(x(t + d))
%                                    = D(d)*f(x(t))
% D*J      - delayed Jacobian        = df(t + d)/dx(t)
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
%__________________________________________________________________________
% %W% Karl Friston %E%

% get dimensions and configure state variables
%--------------------------------------------------------------------------
t     = x(1);                       % peristimulus time (sec)
x     = x(2:end);                   % neuronal states
m     = length(u);  				% number of inputs
n     = length(P.A{1});             % number of sources
x     = reshape(x,n,9);             % neuronal states

% dfdx = [] if t exceeds trial duration (invoking a return to initial state)
%--------------------------------------------------------------------------
if nargout == 1 & (t - P.U) > 1e-6, f = []; return, end

% effective extrinsic connectivity
%--------------------------------------------------------------------------
for i = 1:m
	  P.A{1} = P.A{1} + u(i)*P.B{i};              % forward  connections
      P.A{2} = P.A{2} + u(i)*P.B{i};              % backward connections
      P.A{3} = P.A{3} + u(i)*P.B{i};              % lateral  connections
end

% [default] fixed parameters
%--------------------------------------------------------------------------

E = [32 16 4];           % extrinsic rates (forward, backward, lateral)
G = [1 4/5 1/4 1/4]*128; % intrinsic rates (g1, g2 g3, g4)
D = [2 16];              % propogation delays (intrinsic, extrinsic)
H = [4 32];              % receptor densities (excitatory, inhibitory)
T = [8 16];              % synaptic constants (excitatory, inhibitory)
R = [2 1]/3;             % parameters of static nonlinearity


% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*E(1);
A{2} = exp(P.A{2})*E(2);
A{3} = exp(P.A{3})*E(3);
C    = exp(P.C);

% intrinsic connectivity and parameters
%--------------------------------------------------------------------------
Te   = T(1)/1000*exp(P.T);           % excitatory time constants
Ti   = T(2)/1000;                    % inhibitory time constants
He   = H(1)*exp(P.H);                % excitatory receptor density
Hi   = H(2);                         % inhibitory receptor density

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R    = R.*exp(P.S);
S    = 1./(1 + exp(-R(1)*(x - R(2)))) - 1./(1 + exp(R(1)*R(2)));
dSdx = 1./(1 + exp(-R(1)*(x - R(2)))).^2.*(R(1)*exp(-R(1)*(x - R(2))));

% input
%--------------------------------------------------------------------------
U    = C*spm_erp_u(P,t);

% State: f(x)
%===========================================================================

% Supragranular layer (inhibitory interneurons): Voltage & depolarizing current
%--------------------------------------------------------------------------
f(:,7) = x(:,8);
f(:,8) = (He.*((A{2} + A{3})*S(:,9) + G(3)*S(:,9)) - 2*x(:,8) - x(:,7)./Te)./Te;

% Granular layer (spiny stellate cells): Voltage & depolarizing current
%--------------------------------------------------------------------------
f(:,1) = x(:,4);
f(:,4) = (He.*((A{1} + A{3})*S(:,9) + G(1)*S(:,9) + U) - 2*x(:,4) - x(:,1)./Te)./Te;

% Infra-granular layer (pyramidal cells): depolarizing current
%--------------------------------------------------------------------------
f(:,2) = x(:,5);
f(:,5) = (He.*((A{2} + A{3})*S(:,9) + G(2).*S(:,1)) - 2*x(:,5) - x(:,2)./Te)./Te;

% Infra-granular layer (pyramidal cells): hyperpolarizing current
%--------------------------------------------------------------------------
f(:,3) = x(:,6);
f(:,6) = (Hi*G(4)*S(:,7) - 2*x(:,6) - x(:,3)/Ti)/Ti;

% Infra-granular layer (pyramidal cells): Voltage
%--------------------------------------------------------------------------
f(:,9) = x(:,5) - x(:,6);

% augment with time
%--------------------------------------------------------------------------
f      = [1; f(:)];

if nargout == 1, return, end

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
E  = (A{2} + A{3})*diag(dSdx(:,9)) + G(3)*diag(dSdx(:,9));
E  = diag(He./Te)*E;
S  = sparse(8,9,1,9,9); J = J + kron(S,E);

% Granular layer (spiny stellate cells)
%--------------------------------------------------------------------------
E  = (A{1} + A{3})*diag(dSdx(:,9)) + G(1)*diag(dSdx(:,9));
E  = diag(He./Te)*E;
S  = sparse(4,9,1,9,9); J = J + kron(S,E);

% Infraranular layer (pyramidal cells)
%--------------------------------------------------------------------------
E  = (A{2} + A{3})*diag(dSdx(:,9));
E  = diag(He./Te)*E;
S  = sparse(5,9,1,9,9); J = J + kron(S,E);

E  = G(2)*diag(dSdx(:,1));
E  = diag(He./Te)*E;
S  = sparse(5,1,1,9,9); J = J + kron(S,E);

% Infra-granular layer (pyramidal cells)
%--------------------------------------------------------------------------
E  = G(4)*diag(dSdx(:,7));
E  = Hi*E/Ti;
S  = sparse(6,7,1,9,9); J = J + kron(S,E);

% augment with time
%--------------------------------------------------------------------------
J  = blkdiag(0,J);

% delays
%==========================================================================
De = D(2).*exp(P.D)/1000;
Di = D(1)/1000;
De = (eye(n,n) - 1).*De;
Di = (eye(9,9) - 1)*Di;
De = kron(ones(9,9),De);
Di = kron(Di,eye(n,n));

D  = Di + De;
D  = blkdiag(0,D);

% Implement: dx(t)/dt = f(x(t + d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = D*f = D*J*x(t)
%--------------------------------------------------------------------------
D  = inv(speye(length(J)) - D.*J);

% augment J with delay operator, if D is not requested explicitly
%--------------------------------------------------------------------------
if nargout == 2, J = D*J; end




