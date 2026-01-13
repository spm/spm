function Tab = DEM_FIN(SIM)
% FORMAT Tab = DEM_FIN
%
% Demonstration of COVID-19 modelling using variational Laplace (4 groups)
%__________________________________________________________________________
% This demonstration routine illustrates the complex system modelling of
% financial markets using a generic (Helmholtz Hodge) generative model of
% stochastic chaos under the prior are there exists a pullback attractor.
% This endows the dynamics with a compact functional form; especially when
% using second-order approximations (cf., A generalisation of the Laplace
% approximation).
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% set up and get data
%==========================================================================
cd('C:\Users\karl\Dropbox\Fintech')
rng(0)


%% import data
%--------------------------------------------------------------------------
data = importdata('data_yfinance_daily.csv');

EFT  = {'SPY', 'VEA', 'AGG', 'VNQ', 'DBC', 'BIL'};
name = data.textdata(1,2:end);
iy   =  ismember(name,EFT);
iu   = ~ismember(name,EFT);
date = datenum(data.textdata(2:end,1),'dd/mm/yyyy');
U    = data.data(1:end,iu);
V    = data.data(1:end,iy);

% order assets
%--------------------------------------------------------------------------
i    = [2 3 6 5 4 1];
V    = V(:,i);
EFT  = EFT(i);

% remove null indicator variables
%--------------------------------------------------------------------------
d    = find(min(U) > 0);
U    = U(:,d);
name = name(d);


%% simulation parameters
%==========================================================================
if nargin

    % durations
    %----------------------------------------------------------------------
    N    = SIM.N;               % end point (weeks)
    D    = SIM.D;               % depth of training data (weeks)
    nT   = SIM.nT;              % duration of simulation (in weeks)
    dT   = SIM.dT;              % time between rebalancing (in weeks)

    % specify number of indicator states and assets
    %----------------------------------------------------------------------
    n    = SIM.n;               % number of indicator states
    m    = SIM.m;               % number of assets
    d    = SIM.d;               % order of detrending

else

    % defaults
    %----------------------------------------------------------------------
    N    = 52*7;                   % end point (weeks)
    D    = 512;                 % depth of training data (weeks)
    nT   = 52;                  % duration of simulation (in weeks)
    dT   = 4;                   % time between rebalancing (in weeks)

    % specify number of indicator states and assets
    %----------------------------------------------------------------------
    n    = 3;                   % number of indicator states
    m    = 6;                   % number of assets
    d    = 4;                   % order of detrending

end
dt   = 5;                       % sampling interval (e.g., a week)
N    = size(V,1) - N*dt;        % numer of samples (in days)
U    = U(1:N,:);
V    = V(1:N,:);

% decimation matrix: days to weeks
%--------------------------------------------------------------------------
t    = ceil(N/dt);              % numer of weeks
K    = kron(eye(t,t), ones(1,dt)/dt);
K    = K(:,1:N);                % convolution kernel
K    = times(K,1./sum(K,2));    % returning average
T    = t - nT;                  % start of simulation (weeks)
date = date((1:t)*dt - 1);      % dates

% Exchange-Traded Funds
%--------------------------------------------------------------------------
V    = spm_log(V);              % log value
L    = spm_grad(V);             % rate of log return per day
L    = K*L*dt;                  % rate of log return per week

% eigenreduce
%--------------------------------------------------------------------------
U    = spm_log(U);              % log indicators
U    = K*U;                     % decimate
I    = spm_eigenreduce([U L],1:T,d);

% plot data
%--------------------------------------------------------------------------
spm_figure('GetWin','Data'); clf;

% response and explanatory variables
%--------------------------------------------------------------------------
t    = (T - D):T;
I    = [I(:,1:n)];              % indicator variables
L    = [L(:,1:m)];              % rate of log return
Y    = [I(t,:),L(t,:)];         % training data

subplot(3,1,1)
plot(date,I), title('Log indicator variables','FontSize',14)
datetick('x','mmm-yy'), xlabel('time (mmm-yy')

subplot(3,1,2)
plot(V), title('Log value of ETFs','FontSize',14)
xlabel('time (days)'), spm_axis tight

subplot(3,1,3)
plot(L), hold on
plot([T,T],[-1,1]/4,'-.r')
plot([T,T] - D,[-1,1]/4,'-.k')
title('Log return (per week)','FontSize',14)
xlabel('time (weeks)'), spm_axis tight


% system identification
%==========================================================================

% scaling
%--------------------------------------------------------------------------
scale = std(Y);
X     = rdivide(Y,scale);

% get priors over62 model parameters (polynomial coeficients)
%--------------------------------------------------------------------------
K     = 2;                    % order of polynomial expansion plus one
P.Qp  = 32;                   % variance of parameters (solenoidal)
P.Sp  = 32;                   % variance of parameters (surprisal)
P.Rp  = 32;                   % variance of parameters (expectation)

% prior over precision of states and random fluctuations
%--------------------------------------------------------------------------
W     = var(gradient(X')');
W     = diag(64./W);

% constraints on model parameters (polynomial coefficients)
%--------------------------------------------------------------------------
J      = cell(3,3);
J{1,1} = ones(n,n);
J{2,2} = eye(m,m);
J{2,1} = ones(m,n);
J      = full(spm_cat(J));

% coupling from indicator to response
%--------------------------------------------------------------------------
Q      = cell(3,3);
Q{1,1} = zeros(n,n);
Q{2,2} = zeros(m,m);
Q{2,1} = ones(m,n);
Q      = full(spm_cat(Q));

% get priors 
%--------------------------------------------------------------------------
[pE,pC] = spm_NESS_priors(size(J,1),K,P,W,J,Q);

% scaling
%--------------------------------------------------------------------------
pE.scale = scale(:);
pC.scale = pE.scale*0;


% model inversion with Variational Laplace (fitting flow)
%==========================================================================

% get domain of phase-space and polynomial basis set
%--------------------------------------------------------------------------
x0  = X(1,:)';                 % initial state
M.x = x0;                      % expansion point
M.W = W;                       % precision of random fluctuations
M.K = K;                       % order of polynomial expansion
M.L = K;                       % order of polynomial expansion
M.X = X;                       % legacy points in state-space

% model specification
%--------------------------------------------------------------------------
M.EFT  = EFT(1:m);             % asset names
M.Nmax = 128;                  % maximum number of iterations
M.dt   = dt;                   % sampling interval (days)
M.G    = @spm_NESS_gen_lap;    % generative function
M.pE   = pE;                   % prior expectations (parameters)
M.pC   = pC;                   % prior covariances  (parameters)
M.hE   = log(diag(W));         % prior expectation  (of log-precision)
M.hC   = 1/32;                 % prior covariances  (of log-precision)
M.date = date;

% model inversion with Variational Laplace
%--------------------------------------------------------------------------
B      = spm_ness_U(M);        % get basis functions
F      = gradient(M.X')';      % target flow

% posterior over parameters
%--------------------------------------------------------------------------
[Ep,Cp,Eh] = spm_nlsi_GN(M,B,F);

% posterior precision of random fluctuations
%----------------------------------------------------------------------
W      = diag(exp(Eh));
Ep.W   = W;                    % posterior volatility
M.W    = W;                    % update for flow model

% Bayesian belief updating
%----------------------------------------------------------------------
M.pE   = Ep;                   % and parameters of flow
M.pC   = Cp;                   % and covariance of parameters
M.hE   = Eh;                   % and log-precision


% Bayesian model reduction
%--------------------------------------------------------------------------
% DCM.M  = M;
% DCM.Ep = Ep;
% DCM.Cp = Cp;
% DCM    = spm_dcm_bmr_all(DCM,'All','BMS');
% 
% M.pE   = DCM.M.pE;            % reduced prior expectations (parameters)
% M.pC   = DCM.M.pC;            % reduced prior covariances  (parameters)
% M.P    = DCM.Ep;              % reduced prior covariances  (parameters)
% [Ep,Cp,Eh] = spm_nlsi_GN(M,U,F);
%--------------------------------------------------------------------------


% evaluate Jacobian at initial state
%--------------------------------------------------------------------------
fprintf('\nJacobian at initial conditions\n')
J     = full(spm_diff(@spm_fx_NESS,x0,[],Ep,[],1))
spm_figure('GetWin','CVA: smoothing'); clf;
subplot(2,1,1)
imagesc(J), title('Jacobian (causal coupling)','FontSize',14)
xlabel('latent states'),ylabel('latent states')
axis square

% Bayesian model comparison of enslaving
%==========================================================================

% get reduced priors: first-order effects
%--------------------------------------------------------------------------
J      = cell(3,3);
J{1,1} = ones(n,n);
J{2,2} = eye(m,m);
J{2,1} = ones(m,n);
J      = spm_cat(J);

% coupling from indicator to response
%--------------------------------------------------------------------------
Q      = cell(3,3);
Q{1,1} = zeros(n,n);
Q{2,2} = zeros(m,m);
Q{2,1} = zeros(m,n); %%% remove coupling
Q      = spm_cat(Q);

% get priors 
%--------------------------------------------------------------------------
[rE,rC]  = spm_NESS_priors(size(J,1),K,P,W,J,Q);
rE.scale = scale;
rC.scale = spm_zeros(scale);
F        = spm_log_evidence_reduce(Ep,Cp,pE,pC,rE,rC);
fprintf('\nLog-evidence for 1st order coupling: %.2f\n',F)


% get reduced priors: second-order effects
%--------------------------------------------------------------------------
J      = cell(3,3);
J{1,1} = ones(n,n);
J{2,2} = eye(m,m);
J{2,1} = zeros(m,n); %%% remove coupling
J      = spm_cat(J);

% coupling from indicator to response
%--------------------------------------------------------------------------
Q      = cell(3,3);
Q{1,1} = zeros(n,n);
Q{2,2} = zeros(m,m);
Q{2,1} = ones(m,n);
Q      = spm_cat(Q);

% get priors 
%--------------------------------------------------------------------------
[rE,rC]  = spm_NESS_priors(size(J,1),K,P,W,J,Q);
rE.scale = scale;
rC.scale = spm_zeros(scale);
F        = spm_log_evidence_reduce(Ep,Cp,pE,pC,rE,rC);
fprintf('\nLog-evidence for 2nd order coupling: %.2f\n',F)


% get reduced priors: predicability of indicator variables
%--------------------------------------------------------------------------
J      = cell(3,3);
J{1,1} = eye(n,n); %%% remove coupling
J{2,2} = eye(m,m);
J{2,1} = ones(m,n);
J      = spm_cat(J);

% coupling from indicator to response
%--------------------------------------------------------------------------
Q      = cell(3,3);
Q{1,1} = zeros(n,n);
Q{2,2} = zeros(m,m);
Q{2,1} = ones(m,n);
Q      = spm_cat(Q);

% get priors 
%--------------------------------------------------------------------------
[rE,rC]  = spm_NESS_priors(size(J,1),K,P,W,J,Q);
rE.scale = scale;
rC.scale = spm_zeros(scale);
F        = spm_log_evidence_reduce(Ep,Cp,pE,pC,rE,rC);
fprintf('\nLog-evidence for stochastic chaos: %.2f\n',F)



% model inversion with generlized filtering
%==========================================================================
DEM.EFT = M.EFT;

% serial correlations and orders of motion
%--------------------------------------------------------------------------
DEM.M(1).E.s      = 1/2;       % smoothness of fluctuations
DEM.M(1).E.n      = 1;         % order of generalised motion (states)
DEM.M(1).E.d      = 1;         % order of generalised motion (data)
DEM.M(1).E.nD     = 2;         % number of integrations per dt
DEM.M(1).E.nE     = 2;         % number of iterations
DEM.M(1).E.nN     = 2;         % number of iterations
DEM.M(1).E.linear = 4;         % differentiation scheme

% model
%--------------------------------------------------------------------------
DEM.M(1).f  = @spm_fx_NESS;    % flow
DEM.M(1).g  = @(x,u,P) x.*P.scale; % observer function
DEM.M(1).x  = x0;              % intial state
DEM.M(1).pE = Ep;              % posterior esimates from VL
DEM.M(1).pC = Cp;              % posterior esimates from VL
DEM.M(1).V  = exp(16);
DEM.M(1).W  = W;

% invert
%--------------------------------------------------------------------------
% DEM.Y = Y';
% DEM.U = 1:T;
% DEM   = spm_LAP(DEM);
% 
% % illustrate results
% %------------------------------------------------------------------------
% spm_DEM_qU(DEM.qU)
% subplot(2,2,1), hold on, plot(DEM.Y','.k'), hold off



%% forecasting - analytic
%==========================================================================
spm_figure('GetWin','Fokker Planck'); clf

Nf         = max(2*dT,12);     % forecasting period
DEM.T      = T;                % initial time point
DEM.Y      = Y';               % legacy (training) data
DEM.G      = M;                % generative model of flow
DEM.M(1).x = X(end,:)';        % intial state

% prepare basis functions for analytic forecasting (U)
%------------------------------------------------------------------
N     = numel(x0);
x     = randn(N*N,N);
A.K   = 3;
A.L   = 3;
A.W   = M.W;
DEM.B = spm_ness_U(A,x);

[Ez,Cz,Vz] = spm_NESS_forecast(DEM,Nf);

% overlay actual outcomes if specified
%--------------------------------------------------------------------------
xLim  = [(T - 4*Nf),(T + Nf)];
t     = (1:Nf) + T;
for i = 1:n
    subplot(6,2,i), hold on, set(gca,'ColorOrderIndex',i),
    plot(t,I(t,i),'.')
    set(gca,'XLim',xLim)
end

% enslaved outcomes
%--------------------------------------------------------------------------
for i = (1:m) + n
    subplot(6,2,i), hold on, set(gca,'ColorOrderIndex',i),
    plot(t,L(t,i - n),'.')
    set(gca,'XLim',xLim)
end



%% forecasting - numerical
%==========================================================================
spm_figure('GetWin','Forecasting'); clf

[Ez,Cz,Vz,Py] = spm_NESS_forecasting(DEM,Nf,n,m);

% overlay actual outcomes if specified
%--------------------------------------------------------------------------
t  = (1:Nf) + T;
subplot(2,2,1), hold on, set(gca,'ColorOrderIndex',1),
for i = 1:n
    plot(t,I(t,i),'.')
end

% overlay actual outcomes if specified
%--------------------------------------------------------------------------
xLim  = [(T - 4*Nf),(T + Nf)];
t     = (1:Nf) + T;
for i = 1:n
    subplot(8,3,12 + i), hold on, set(gca,'ColorOrderIndex',i),
    plot(t,I(t,i),'.')
    set(gca,'XLim',xLim)
end

% enslaved outcomes
%--------------------------------------------------------------------------
for i = (1:m) + n
    subplot(8,3,12 + i), hold on, set(gca,'ColorOrderIndex',i),
    plot(t,L(t,i - n),'.')
    set(gca,'XLim',xLim)
end

if ~m, return, end

%% Projections
%==========================================================================
spm_figure('GetWin','Responses'); clf

% extract response variables
%--------------------------------------------------------------------------
z     = (1:nT) + T;            % future time points with data
t     = (1:Nf) + T - 1;        % future time points predicted
i     = (1:m) + n;             % indices of response variables
EL    = Ez(i,:);               % posterior predictive expectations
VL    = Vz(i,:);               % posterior predictive variances

% plot credible intervals and data
%--------------------------------------------------------------------------
xLim  = [T - 2*Nf,T + Nf];
for i = 1:m

    % predictive posterior
    %----------------------------------------------------------------------
    subplot(6,2,i), hold off
    set(gca,'ColorOrderIndex',i + n)
    spm_plot_ci(EL(i,:),VL(i,:),t), hold on

    % overlay past and future outcomes
    %----------------------------------------------------------------------
    set(gca,'ColorOrderIndex',i + n)
    plot(1:T,L(1:T,i),'LineWidth',2), hold on
    set(gca,'ColorOrderIndex',i + n)
    plot(z,L(z,i),'.'), plot(xLim,[0,0],':'),hold off
    title(['RoR: ' EFT{i}],'FontSize',14)
    set(gca,'XLim',xLim)

end

% simulate portfolio management under this genertive model
%==========================================================================
spm_figure('GetWin','Portfolio forecasts'); clf

% constraints on Bayesian belief updating
%--------------------------------------------------------------------------
Rp         = pC;
Rp.Qp      = spm_zeros(pC.Qp);
Rp.Sp      = spm_zeros(pC.Sp);

Rp.Rp      = (pC.Rp > 0)/6;
Rp         = diag(spm_vec(Rp));
DEM.G.pC   = Rp*DEM.M.pC*Rp;

DEM.RoR    = 60;
DEM.Sharpe = 1;
Tab        = spm_portfolio(L,I,DEM,T,dT)

return

% % line search for Rp.Rp
% %--------------------------------------------------------------------------
% for i = 1:8
%     Rp     = pC;
%     Rp.Qp  = spm_zeros(pC.Qp);
%     Rp.Sp  = spm_zeros(pC.Sp);
%     Rp.Rp  = (pC.Rp > 0)/i;
%     Rp     = diag(spm_vec(Rp));
%     M.pC   = Rp*DEM.M.pC*Rp;
% 
%     spm_figure('GetWin','Portfolio forecasts'); clf
%     Tab    = spm_portfolio(L,I,M,DEM,T,dT)
%     tab    = table2array(Tab);
%     RoR(i) = tab(5,1)
%     Rat(i) = tab(5,3)
% end

% line search for expected utility (RoR) under a Sharpe ratio of 2
%--------------------------------------------------------------------------
% for i = 1:8
%     DEM.RoR    = 20*i;
%     DEM.Sharpe = 1/2;
%     spm_figure('GetWin','Portfolio forecasts'); clf
%     Tab    = spm_portfolio(L,I,M,DEM,T,dT)
%     tab    = table2array(Tab);
%     RoR(i) = tab(5,1)
%     Rat(i) = tab(5,3)
% end

% repeat for increasing intervals between rebalancing
%==========================================================================
% for i = 1:11
%     spm_figure('GetWin','Portfolio forecasts'); clf
%     Tab    = spm_portfolio(L,I,M,DEM,T,i + 1)
%     tab    = table2array(Tab);
%     RoR(i) = tab(5,1);
%     Rat(i) = tab(5,3);
% end
% 
% spm_figure('GetWin','Re-balancing'); clf
% subplot(2,2,1), bar((1:numel(RoR)) + 1,RoR)
% title('Annualised rate of return','FontSize',14)
% xlabel('number of weeks'), ylabel('%')
% axis square, box off
% 
% subplot(2,2,2), bar((1:numel(Rat)) + 1,Rat)
% title('Sharpe ratio','FontSize',14)
% xlabel('number of weeks'), ylabel('ratio')
% axis square, box off
%  
% return

% free energy and Market cycles
%==========================================================================
spm_figure('GetWin','Financial free energy'); clf

% surprisal in the past
%--------------------------------------------------------------------------
F     = [];
d     = 128;
for t = (T - d):(T + dT - 1)
    x     = [I(t,:) L(t,:)]./scale;
    [~,S] = spm_NESS_gen_lap(Ep,M,x);
    F(t)  = S;
end

% surprisal forecasts
%--------------------------------------------------------------------------
Fi    = [];
for i = 1:size(Py,3)
    f = [];
    for t = 1:dT
        x     = Py(:,t,i)'./scale;
        [~,s] = spm_NESS_gen_lap(Ep,M,x);
        f(t)  = s;
    end

    % prepend legacy data
    %----------------------------------------------------------------------
    Fi(i,:) = [F(1:T) f(2:end)];
end

% free energy fluctuations
%--------------------------------------------------------------------------
subplot(2,1,1), hold off
t = M.date;
i = (T - d):(T + dT - 1);
plot(t(i),Fi(:,i)','r','LineWidth',1/4), hold on
plot(t(i),F(i),'b','LineWidth',1)
datetick('x','mmm-yy'), xlabel('time'), ylabel('surprisal (nats)');
title('Financial free energy','FontSize',14)

% Phase portraits
%--------------------------------------------------------------------------
S     = interp(F,32);
for i = 1:size(Fi,1)

    % interpolate free energy fluctuations
    %----------------------------------------------------------------------
    Si(i,:) = interp(Fi(i,:),32);
end

subplot(2,1,2), hold off
dSdt  = gradient(S);
dFdt  = gradient(Si);
plot(dFdt',Si','r:','LineWidth',1/4), hold on
plot(dSdt,S,'b','LineWidth',1)
plot([0,0],get(gca,'YLim'),':k')
axis square, title('Phase portrait','FontSize',14)
xlabel('time derivative'), ylabel('surprisal (nats)');

return

% subroutines
%==========================================================================

function X = spm_grad(X)
% gradients based on diff operator (i.e. the past)
% X  - mattrix of variables
%__________________________________________________________________________

X  = diff(X);
X  = [X(1,:); X];

return


function I = spm_eigenreduce(U,T,d)
% eigenreduction of indicator variables
% FORMAT I = spm_eigenreduce(U,T,[d])
%--------------------------------------------------------------------------
% U  - indicator variables
% T  - based on first T samples
% d  - order of DCT (detrending)
%__________________________________________________________________________
% U  = spm_detrend(U);
% U  = spm_en(U);
% U  = spm_en(U')';

n  = size(U,1);
X  = spm_dctmtx(n,d);
B  = X(T,:)\U(T,:);
U  = U - X*B;

% get eigenvectors
%--------------------------------------------------------------------------
v  = spm_svd(spm_en(U(T,:))');
I  = U*v;
s  = diag(1./std(I(T,:)));
I  = I*s;

return

function [X] = spm_exp_conv(X,t)
% Epponential convolution
% FORMAT [X] = spm_exp_conv(X,t)
% X    - matrix
% sx   - kernel length
%__________________________________________________________________________
%
% spm_exp_conv is a one dimensional convolution of a matrix variable in
% working memory.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1999-2022 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
n  = size(X,1);
K  = exp(-(0:(6*t)).^2/(2*(t^2)));
C  = spm_convmtx(K(:),n);
C  = C(1:n,1:n);
C  = diag(1./sum(C,2))*C;
X  = C*X;

return


% Iterated simulations
%==========================================================================
tab   = zeros(5,4,8);
for i = 1:8

    % durations
    %----------------------------------------------------------------------
    SIM.N  = (i - 1)*52;        % end point (weeks)
    SIM.D  = 512;               % depth of training data (weeks)
    SIM.nT = 52;                % duration of simulation (in weeks)
    SIM.dT = 4;                 % time between rebalancing (in weeks)

    % specify number of indicator states and assets
    %----------------------------------------------------------------------
    SIM.n  = 3;                 % number of indicator states
    SIM.m  = 6;                 % number of assets
    SIM.d  = 4;                 % order of detrending

    Tab = DEM_FIN(SIM);
    tab(:,:,i) = table2array(Tab);

end

% bar chart results
%--------------------------------------------------------------------------
spm_figure('GetWin','Annual performance'); clf
for i = 1:4
    subplot(2,2,i)
    bar(squeeze(tab(:,i,:))')
    title(Tab.Properties.VariableNames{i})
    legend(Tab.Properties.RowNames)
    axis square
end


% NOTES: numerical checks on Jacobian
%==========================================================================

% specify
%--------------------------------------------------------------------------
n    = 3;                       % number of enslaving states
m    = 2;                       % number of enslaved  states
K    = 2;                       % order of polynomial expansion
W    = eye(n + m,n + m);

% constraints on model parameters (polynomial coefficients)
%--------------------------------------------------------------------------
J      = cell(3,3);
J{1,1} = ones(n,n);
J{2,2} = eye(m,m);
J{2,1} = ones(m,n);
J      = full(spm_cat(J));

% coupling from indicator to response
%--------------------------------------------------------------------------
Q      = cell(3,3);
Q{1,1} = zeros(n,n);
Q{2,2} = zeros(m,m);
Q{2,1} = ones(m,n);
Q      = full(spm_cat(Q));

% get priors 
%--------------------------------------------------------------------------
[pE,pC] = spm_NESS_priors(size(J,1),K,1,W,J,Q);

% numerical check of Jacobian
%--------------------------------------------------------------------------
Ep   = spm_unvec(spm_vec(pC)/2,pE);
Ep.W = W;
M.W  = W;
M.K  = K;
M.L  = K;
x    = ones(n + m,1);
full(spm_diff(@spm_fx_NESS,x,[],Ep,[],1))

% ensure numerical consistency of flow
%--------------------------------------------------------------------------
spm_NESS_gen_lap(Ep,M,x)'
spm_fx_NESS(x,[],Ep)

% ensure numerical equivalence and consistency of surprisal gradients
%--------------------------------------------------------------------------
spm_fx_NESS(x,[],Ep,[],'DS')
full(spm_diff(@spm_fx_NESS,x,[],Ep,[],'S',1)')

cell2mat(spm_NESS_gen_lap(Ep,M,x,'DS'))
full(spm_diff(@spm_NESS_gen_lap,Ep,M,x,'S',3)')

% effect of parameters on Jacobian
%--------------------------------------------------------------------------
dJdP  = spm_diff(@spm_fx_NESS,x,[],Ep,[],[1 3]);
