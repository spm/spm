function DEM_ness_Lap
% Nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT DEM_ness_Lap
%--------------------------------------------------------------------------
% This demonstration routine illustrates the use of generalised
% (variational) filtering to invert a generic model of any complex
% dynamical system. This model makes the minimal assumption that there
% exists a nonequilibrium steady-state (i.e., a pullback attractor). The
% ensuing solution to the density dynamics—under the Helmholtz-Hodge
% decomposition—furnishes a normal form for the flow or dynamics. This
% normal form can be parameterised to 2nd order in terms of (solenoidal and
% dissipative) flow operators that operate on the gradients of the system’s
% potential, surprisal or self-information (i.e., the implausibility of
% being found in any given state and the NESS density). The second order
% parameterisation (of the underlying flow) is used as a generative model
% to learn the parameters of the system.
%
% This enables forecasting; in the sense that the system can be solved for
% future time points (with random fluctuations on the flow) to produce an
% ensemble or sample estimate of future states.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%% Create a model under polynomial NESS
%==========================================================================
rng(0)

% states and precisons
%--------------------------------------------------------------------------
nm  = 3;                      % number of master states
ns  = 1;                      % number of slave states
n   = nm + ns;                % total number of states
K   = 2;                      % order of polynomial expansion (suprisal)
L   = 2;                      % order of polynomial expansion (solenoidal)
V   = 32;                     % precision of parameters (prior)
W   = 256;                    % precision of random fluctuations

% times
%--------------------------------------------------------------------------
dt   = 1;                     % sampling interval (e.g., days)
nT   = 32;                    % duration of forecasting period
T    = 128;                   % number of past samples

% constraints on model parameters (polynomial coeficients)
%--------------------------------------------------------------------------
i      = (1:ns) + nm;
j      = 1:nm;
J      = ones(n,n);           % contraints on Jacobian
J(j,i) = 0;                   % i.e., causal coupling
J(i,j) = 0;                   % i.e., causal coupling

R      = zeros(n,n);          % enslaving state (first master state)
R(i,1) = 1;                   % and enslaved states

W      = ones(1,n)*W;         % precision of random fluctuations
W(i)   = 256;                 % for enslaved state
W      = diag(W);

C      = ones(1,n);            % covaraince of NESS density
C(i)   = 1/8;                  % for enslaved state
C      = diag(C);

% get priors and sample some parameters (i.e., polynomial coefficients)
%--------------------------------------------------------------------------
[pE,pC] = spm_NESS_priors(n,K,V,W,J,R,C);
pC      = diag(spm_vec(pC));
Np      = spm_length(pE);
P       = spm_unvec(spm_vec(pE) + sqrt(pC)*randn(Np,1)/4,pE);
P.Rp    = ~~P.Rp/2;
pP.P{1} = P;

% evaluate Jacobian at initial state
%--------------------------------------------------------------------------
disp(full(spm_diff(@spm_fx_NESS,ones(n,1),[],P,[],1)))


%% solution of stochastic and underlying ordinary differential equation
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS'); clf

x0  = randn(n,1);             % initial state
M.f = @spm_fx_NESS;           % flow under NESS
M.x = x0;                     % expansion point

% deterministic solution
%--------------------------------------------------------------------------
U.dt = dt;
U.u  = zeros(T + nT,1);
Z    = spm_int_ode(P,M,U);
Y.y  = Z(1:T,:);

subplot(3,1,1)
plot(Z), hold on, set(gca,'ColorOrderIndex',1)
title('Trajectory: deterministic','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off


% stochastic solution
%--------------------------------------------------------------------------
DEM.M(1).f  = @spm_fx_NESS;    % dynamics
DEM.M(1).g  = @(x,u,P) x;      % observer function
DEM.M(1).x  = x0;              % inital state
DEM.M(1).pE = P;               % parameters
DEM.M(1).V  = exp(8);          % precision of data
DEM.M(1).W  = W;               % precision of dynamics
 
% orders of generalised motion
%--------------------------------------------------------------------------
DEM.M(1).E.s      = 1/2;       % smoothness of fluctuations
DEM.M(1).E.n      = 2;         % order of generalised motion (states)
DEM.M(1).E.d      = 1;         % order of generalised motion (data)
DEM.M(1).E.nD     = 2;         % number of integrations per dt


DEM   = spm_DEM_generate(DEM.M,T + nT);
Z     = DEM.Y';                % past and future outcomes
Y.y   = Z(1:T,:);              % past data

subplot(3,1,2)
plot(Z), hold on, set(gca,'ColorOrderIndex',1)
title('Trajectory: stochastic','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off

% NESS density and expected flow
%--------------------------------------------------------------------------
M.W = W;
M.K = K;
M.L = L;

[F,S,Q,L,H] = spm_NESS_gen_lap(P,M,x0);


% model inversion with Variational Laplace (fitting motion)
%==========================================================================

% get domain of phase-space and polynomial basis set
%--------------------------------------------------------------------------
M      = rmfield(M,'f');
M.X    = Y.y;
M.W    = W;
U      = spm_ness_U(M);        % get state space and flow
F      = gradient(M.X',dt)';   % target flow

% model specification
%--------------------------------------------------------------------------
M.Nmax = 64;                   % maximum number of iterations
M.G    = @spm_NESS_gen_lap;    % generative function
M.pE   = pE;                   % prior expectations (parameters)
M.pC   = pC;                   % prior covariances  (parameters)
M.hE   = log(diag(W));         % prior expectation  (of log-precision)
M.hC   = 1/32;                 % prior covariances  (of log-precision)

% model inversion with Variational Laplace
%--------------------------------------------------------------------------
[Ep,Cp,Eh] = spm_nlsi_GN(M,U,F);
W          = diag(exp(Eh));
Ep.W       = W;

% model inversion with generlized filtering
%==========================================================================

% serial correlations and orders of motion
%--------------------------------------------------------------------------
DEM.M(1).E.s      = 1/2;       % smoothness of fluctuations
DEM.M(1).E.n      = 2;         % order of generalised motion (states)
DEM.M(1).E.d      = 1;         % order of generalised motion (data)
DEM.M(1).E.nD     = 2;         % number of integrations per dt
DEM.M(1).E.nN     = 2;         % number of interations
DEM.M(1).E.nE     = 8;         % number of interations
DEM.M(1).E.linear = 1;         % differentiation scheme

% model
%--------------------------------------------------------------------------
DEM.M(1).f  = @spm_fx_NESS;    % dynamics
DEM.M(1).g  = @(x,u,P) x;      % observer function
DEM.M(1).x  = x0;              % intial state
DEM.M(1).pE = Ep;              % posterior esimates from VL
DEM.M(1).pC = Cp;              % posterior esimates from VL
DEM.M(1).V  = exp(16);
DEM.M(1).W  = W;

% invert
%--------------------------------------------------------------------------
DEM.Y = Y.y';
DEM.U = zeros(1,T);
DEM   = spm_DEM(DEM);

% illustrate results
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU)
subplot(2,2,1), hold on, plot(DEM.Y','.k'), hold off

spm_figure('GetWin','parameters'); clf
spm_DEM_qP(DEM.qP,pP)



%% forecasting
%==========================================================================
spm_figure('GetWin','forecasting'); clf

spm_NESS_forecasting(DEM,nT,nm,ns);

% overlay actual outcomes if specified
%--------------------------------------------------------------------------
z = (1:nT) + T;
subplot(2,2,1), hold on, set(gca,'ColorOrderIndex',1),
for i = 1:size(Z,2)
    plot(z,Z(z,i),'.')
end

% enslaving outcomes
%--------------------------------------------------------------------------
subplot(4,1,3), hold on, set(gca,'ColorOrderIndex',1),
for i = 1:1:nm
    plot(z,Z(z,i),'.')
end

% enslaved outcomes
%--------------------------------------------------------------------------
for i = (nm + 1):n
    subplot(4,1,4), hold on
    set(gca,'ColorOrderIndex',i),
    plot(z,Z(z,i),'.')
end

return



% NOTES

function NESS = spm_Lorenz2Lap
% FORMAT NESS = spm_Lorenz2Lap
% return a numerical estimate of polynomial coefficients for Lorenz system
%__________________________________________________________________________

%% Illustration of high order density learning
%==========================================================================
% dxdt = f(x) + w:  see notes at the end of this script
%--------------------------------------------------------------------------
f    = @(x,v,P,M) [P(1)*x(2) - P(1)*x(1);
                   P(3)*x(1) - x(2) - x(1)*x(3);
                  -P(2)*x(3) + x(1)*x(2)]/64;
J    = @(x,v,P,M) [[     -P(1),P(1),     0];
                   [P(3) - x(3), -1, -x(1)];
                   [     x(2), x(1), -P(2)]]/64;
P    = [10; 8/3; 32];                  % parameters [sig, beta, rho]
x0   = [1; 1; 24];                     % initial state

% state-space model (for SPM integrators)
%--------------------------------------------------------------------------
M.f  = f;
M.J  = J;
M.g  = @(x,v,P,M) x;
M.x  = x0;
M.m  = 0;
M.pE = P;
M.W  = diag([1/8 1/16 1/32]);           % precision of random fluctuations
M.K  = 2;                               % order of expansion (suprisal)
M.L  = 3;                               % order of expansion (flow)

% state-space for (Laplace) solution 
%--------------------------------------------------------------------------
N    = 8;                                % number of bins
d    = 12;                               % distance
m    = [0 0 28];
x{1} = linspace(m(1) - d,m(1) + d,N);
x{2} = linspace(m(2) - d,m(2) + d,N);
x{3} = linspace(m(3) - d,m(3) + d,N);

% estimate
%--------------------------------------------------------------------------
NESS = spm_ness_Lap(M,x);



