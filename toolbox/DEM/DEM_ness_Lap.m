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
rng(1)

% states and precisons
%--------------------------------------------------------------------------
n   = 3;                      % number of states
K   = 2;                      % order of polynomial expansion (suprisal)
V   = 128;                     % precision of parameters (prior)
W   = 256;                    % precision of random fluctuations

% times
%--------------------------------------------------------------------------
dt   = 1;                     % sampling interval (e.g., days)
nT   = 64;                    % duration of forecasting period
T    = 128;                   % number of past samples

% get model parameters (polynomial coeficients)
%--------------------------------------------------------------------------
[pE,pC] = spm_NESS_priors(n,K,V,W);
Np      = spm_length(pE);
P       = spm_unvec(spm_vec(pE) + sqrt(pC)*randn(Np,1)/4,pE);
pP.P{1} = P;             


%% solution of stochastic and underlying ordinary differential equation
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS'); clf

x0  = randn(n,1);             % initial state
M.f = @spm_fx_NESS;           % flow under NESS
M.x = x0;                     % expansion point
M.W = W;                      % precision of random fluctuations
M.K = K;                      % order of polynomial expansion (suprisal)

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
DEM.M(1).V  = 512;             % precision of data
DEM.M(1).W  = W;               % precision of dynamics
 
% orders of generalised motion
%--------------------------------------------------------------------------
DEM.M(1).E.s      = 1/2;       % smoothness of fluctuations
DEM.M(1).E.n      = 2;         % order of generalised motion (states)
DEM.M(1).E.d      = 1;         % order of generalised motion (data)
DEM.M(1).E.nD     = 2;         % number of integrations per dt
DEM.M(1).E.linear = 4;         % differentiation scheme


DEM   = spm_DEM_generate(DEM.M,T + nT);
Z     = DEM.Y';                % past and future outcomes
Y.y   = Z(1:T,:);              % past data

subplot(3,1,2)
plot(Z), hold on, set(gca,'ColorOrderIndex',1)
title('Trajectory: stochastic','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off

%%% return

% model inversion with Variational Laplace (fitting motion)
%==========================================================================

% get domain of phase-space and polynomial basis set
%--------------------------------------------------------------------------
M      = rmfield(M,'f');
M.X    = Y.y;
U      = spm_ness_U(M);        % get state space and flow
F      = gradient(M.X',dt)';   % target flow

% model specification
%--------------------------------------------------------------------------
M.Nmax = 64;                   % maximum number of iterations
M.G    = @spm_NESS_gen_lap;    % generative function
M.pE   = pE;                   % prior expectations (parameters)
M.pC   = pC;                   % prior covariances  (parameters)
M.hE   = log(W);               % prior expectation  (of log-precision)
M.hC   = 1/512;                % prior covariances  (of log-precision)

% model inversion with Variational Laplace
%--------------------------------------------------------------------------
[Ep,Cp] = spm_nlsi_GN(M,U,F);

% NESS density and expected flow
%--------------------------------------------------------------------------
% [F,S,Q,L,H] = spm_NESS_gen_lap(Ep,M,U);


% model inversion with generlized filtering
%==========================================================================

% serial correlations and orders of motion
%--------------------------------------------------------------------------
DEM.M(1).E.s      = 1/2;       % smoothness of fluctuations
DEM.M(1).E.n      = 2;         % order of generalised motion (states)
DEM.M(1).E.d      = 1;         % order of generalised motion (data)
DEM.M(1).E.nD     = 2;         % number of integrations per dt
DEM.M(1).E.linear = 4;         % differentiation scheme

% model
%--------------------------------------------------------------------------
DEM.M(1).f  = @spm_fx_NESS;    % dynamics
DEM.M(1).g  = @(x,u,P) x;      % observer function
DEM.M(1).x  = x0;              % intial state
DEM.M(1).pE = Ep;              % posterior esimates from VL
DEM.M(1).pC = pC;              % posterior esimates from VL
DEM.M(1).V  = 512;
DEM.M(1).W  = W;

% invert
%--------------------------------------------------------------------------
DEM.Y = Y.y';
DEM.U = zeros(1,T);
DEM   = spm_LAP(DEM);

% illustrate results
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU)
subplot(2,2,1), hold on, plot(DEM.Y','.k'), hold off

spm_figure('GetWin','parameters'); clf
spm_DEM_qP(DEM.qP,pP)



%% forecasting
%==========================================================================
spm_figure('GetWin','forecasting'); clf

spm_NESS_forecasting(DEM,nT);

% overlay actual outcomes if specified
%--------------------------------------------------------------------------
z = (1:nT) + T;
subplot(2,2,1), hold on
plot(z,Z(z,:),'.k'), hold off

subplot(2,1,2), hold on
plot(z,Z(z,:),'.k'), hold off

return




