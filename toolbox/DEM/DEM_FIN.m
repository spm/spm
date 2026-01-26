function [Tab,F,DEM] = DEM_FIN(SIM,Ep,ARG)
% FORMAT [Tab,F,DEM] = DEM_FIN(SIM,ARG)
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
V    = data.data(1:end,iy);     % values of assets (ETF)
O    = data.data(1:end,iu);     % other (indicator) variables

% order assets (risk-free first)
%--------------------------------------------------------------------------
i    = [2 3 6 5 4 1];
V    = V(:,i);
EFT  = EFT(i);

% remove null indicator variables
%--------------------------------------------------------------------------
d    = find(min(O) > 0);
O    = O(:,d);
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
    N    = 52*(1 - 1);          % end point (weeks)
    D    = 256;                 % depth of training data (weeks)
    nT   = 52;                  % duration of simulation (in weeks)
    dT   = 4;                   % time between rebalancing (in weeks)

    % specify number of indicator states and assets
    %----------------------------------------------------------------------
    n    = 3;                   % number of indicator states
    m    = 6;                   % number of assets
    d    = 8;                   % order of detrending

end
dt   = 5;                       % sampling interval (e.g., a week)
N    = size(V,1) - N*dt;        % number of samples (in days)
O    = O(1:N,:);                % indicator variables
V    = V(1:N,:);                % value variables


% decimation matrix: days to weeks
%--------------------------------------------------------------------------
t    = ceil(N/dt);              % number of weeks
K    = kron(eye(t,t), ones(1,dt)/dt);
K    = K(:,1:N);                % convolution kernel
K    = times(K,1./sum(K,2));    % returning average
T    = t - nT;                  % start of simulation (weeks)
date = date((1:t)*dt - 1);      % dates

% Exchange-Traded Funds
%--------------------------------------------------------------------------
U    = spm_dctmtx(D,d,1:t);     % exogenous inputs
V    = spm_log(V);              % log value
L    = spm_grad(V);             % rate of log return per day
L    = K*L*dt;                  % rate of log return per week

% eigenreduce
%--------------------------------------------------------------------------
O    = spm_log(O);              % log indicators
O    = K*O;                     % decimate
I    = spm_eigenreduce([O L],1:T,3);

% plot data
%--------------------------------------------------------------------------
spm_figure('GetWin','Data'); clf;

% response and explanatory variables
%--------------------------------------------------------------------------
t    = (T - D):T;
I    = [I(:,1:n)];              % indicator variables
L    = [L(:,1:m)];              % rate of log return
Y    = [I(t,:),L(t,:)];         % training data
u    = [U(t,:)];                % exogenous inputs 

subplot(3,1,1)
plot(date,I(:,1:n)), title('Log indicator variables','FontSize',14)
datetick('x','mmm-yy'), xlabel('time (mmm-yy')

subplot(3,1,2)
plot(V), title('Log value of ETFs','FontSize',14)
xlabel('time (days)'), spm_axis tight

subplot(3,1,3)
plot(U,'--'), hold on
plot(L)
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

% get priors over model parameters (polynomial coeficients)
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

% mean
%--------------------------------------------------------------------------
pE.Rp(1,:) = 0;
pC.Rp(1,:) = 0;

% scaling
%--------------------------------------------------------------------------
pE.scale = scale(:);
pC.scale = pE.scale*0;

% scaling
%--------------------------------------------------------------------------
nx   = n + m;
pE.U = zeros(d,nx);
pC.U = zeros(d,nx);
pC.U(:,1:nx) = 32;

% model inversion with Variational Laplace (fitting flow)
%==========================================================================

% get domain of phase-space and polynomial basis set
%--------------------------------------------------------------------------
x0  = X(1,:)';                 % initial state
u0  = u(1,:)';                 % initial cause
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
M.IS   = @spm_ness_F;          % generative function
M.pE   = pE;                   % prior expectations (parameters)
M.pC   = pC;                   % prior covariances  (parameters)
M.hE   = log(diag(W));         % prior expectation  (of log-precision)
M.hC   = 1/32;                 % prior covariances  (of log-precision)
M.date = date;

% model inversion with Variational Laplace
%--------------------------------------------------------------------------
B      = spm_ness_U(M);        % get basis functions
B.u    = u(1:(D + 1),:);       % add exogenous inputs
F      = gradient(M.X')';      % target flow

% posterior over parameters
%--------------------------------------------------------------------------
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,B,F);
F0     = F;

% posterior precision of random fluctuations
%--------------------------------------------------------------------------
W      = diag(exp(Eh));
Ep.W   = W;                    % posterior volatility
M.W    = W;                    % update for flow model

% Bayesian belief updating
%--------------------------------------------------------------------------
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
% [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,F);
%--------------------------------------------------------------------------


% evaluate Jacobian at initial state
%--------------------------------------------------------------------------
spm_figure('GetWin','Jacobian'); clf;
spm_J = @(x,u,P) full(spm_diff(@spm_ness_f,x,u,P,[],1));

fprintf('\nJacobian at initial conditions\n')
J     = spm_J(x0,u0,Ep)
subplot(2,2,1)
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
J{2,1} = zeros(m,n); %%% remove coupling
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
rE.scale = pE.scale;
rC.scale = pC.scale;
rE.U     = pE.U;
rC.U     = pC.U*0;   %%% remove coupling
F1       = spm_log_evidence_reduce(Ep,Cp,pE,pC,rE,rC);
fprintf('\nLog-evidence for enslaving: %.2f\n',F1)

% get reduced priors: predicability of indicator variables
%--------------------------------------------------------------------------
J      = cell(3,3);
J{1,1} = eye(n,n);   %%% remove coupling
J{2,2} = eye(m,m);
J{2,1} = eye(m,n);
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
rE.scale = pE.scale;
rC.scale = pC.scale;
rE.U     = pE.U;
rC.U     = pC.U;
F2       = spm_log_evidence_reduce(Ep,Cp,pE,pC,rE,rC);
fprintf('\nLog-evidence for stochastic chaos: %.2f\n',F2)



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
DEM.M(1).f  = @spm_ness_f;     % flow
DEM.M(1).g  = @(x,u,P) x.*P.scale; % observer function
DEM.M(1).x  = x0;              % intial state
DEM.M(1).pE = Ep;              % posterior esimates from VL
DEM.M(1).pC = Cp;              % posterior esimates from VL
DEM.M(1).V  = exp(16);
DEM.M(1).W  = W;

DEM.M(2).x  = [];
DEM.M(2).v  = u0;              % intial cause
DEM.M(2).V  = exp(16);         % prior precicion


% invert
%--------------------------------------------------------------------------
DEM.Y = Y';
DEM.U = B.u';

% DEM   = spm_LAP(DEM);
% 
% % illustrate results
% %------------------------------------------------------------------------
% spm_DEM_qU(DEM.qU)
% subplot(2,2,1), hold on, plot(DEM.Y','.k'), hold off



%% forecasting - analytic
%==========================================================================
spm_figure('GetWin','Fokker Planck'); clf

Nf         = 16;               % forecasting period
DEM.T      = T;                % initial time point
DEM.Y      = Y';               % legacy (training) data
DEM.G      = M;                % generative model of flow
DEM.M(1).x = X(end,:)';        % intial state
DEM.M(2).v = u(end,:)';        % intial cause

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

% future causes
%--------------------------------------------------------------------------
t     = (1:Nf) + T;
DEM.U = U(t,:)';

% predictive posterior realizations
%--------------------------------------------------------------------------
[Ez,Cz,Vz,Py] = spm_NESS_forecasting(DEM);

% overlay actual outcomes if specified
%--------------------------------------------------------------------------
subplot(2,2,1), hold on, set(gca,'ColorOrderIndex',1),
for i = 1:n
    plot(t,I(t,i),'.')
end

% overlay actual outcomes if specified
%--------------------------------------------------------------------------
xLim  = [(T - 4*Nf),(T + Nf)];
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

% simulate portfolio management under this generative model
%==========================================================================

% constraints on Bayesian belief updating
%--------------------------------------------------------------------------
Rp       = spm_zeros(pC);
Rp.Rp    = (pC.Rp > 0)/8;
Rp       = diag(spm_vec(Rp));
DEM.G.pC = Rp*Cp*Rp;

% simulate portfolio management
%--------------------------------------------------------------------------
Tab = spm_portfolio(L,I,U,DEM,dT)

% outputs if requested
%--------------------------------------------------------------------------
F = [F0,F1,F2];
if nargin > 1
    F = eval(ARG);
end

return

% % line search for Rp.Rp
% %------------------------------------------------------------------------
% for i = 1:16
% 
%     Rp       = pC;
%     Rp.Qp    = spm_zeros(pC.Qp);
%     Rp.Sp    = spm_zeros(pC.Sp);
% 
%     Rp.Rp    = (pC.Rp > 0)/i;
%     Rp.U     = (pC.U  > 0)/i;
%     Rp       = diag(spm_vec(Rp));
%     DEM.G.pC = Rp*Cp*Rp;
% 
%     spm_figure('GetWin','Portfolio forecasts'); clf
%     Tab    = spm_portfolio(L,I,U,DEM,dT)
%     tab    = table2array(Tab);
%     RoR(i) = tab(5,1)
%     Rat(i) = tab(5,3)
% 
% end


% free energy and Market cycles
%==========================================================================
spm_figure('GetWin','Financial free energy'); clf

% surprisal in the past
%--------------------------------------------------------------------------
F     = [];
d     = 128;
for t = (T - d):(T + dT - 1)
    x     = [I(t,:) L(t,:)]./scale;
    u     = U(t,:);
    F(t)  = DEM.M(1).f(x,u,Ep,[],'S');
end

% surprisal forecasts
%--------------------------------------------------------------------------
Fi    = [];
for i = 1:size(Py,3)
    f = [];
    for t = 1:dT
        x     = Py(:,t,i)'./scale;
        u     = U(T + t,:);
        f(t)  = DEM.M(1).f(x,u,Ep,[],'S');
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


%% Is the Market chaotic?
%==========================================================================
spm_figure('GetWin','Jacobian'); clf

% state space for evaluation
%--------------------------------------------------------------------------
ni    = 64;
nx    = n + m;
E     = zeros(nx,ni);
F     = zeros(nx,ni);
J     = zeros(nx,nx,ni);
for i = 1:ni
    F(:,i)   = spm_ness_f(X(i,:),u(i,:),Ep);
    J(:,:,i) = spm_J(X(i,:),u(i,:),Ep);
    E(:,i)   = sort(real(eig(J(:,:,i))),'descend');
end

subplot(2,2,1)
imagesc(mean(J,3)), title('Jacobian (causal coupling)','FontSize',14)
xlabel('latent states'),ylabel('latent states'), axis square

% Lyapunov exponent and Hausdorff dimension (Kaplan-Yorke conjecture)
%--------------------------------------------------------------------------
LE    = mean(E,2);
j     = sum(LE >= 0);
CD    = j + sum(LE(1:j))/abs(LE(j + 1))

subplot(2,1,2)
i     = 1:ni;
x     = X';
quiver3(x(1,i),x(2,i),x(3,i),F(1,i),F(2,i),F(3,i))
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
title('Flow','Fontsize',16), axis square

subplot(2,2,2), hold off
plot3(x(1,i),x(2,i),x(3,i),'.b'), hold on
i     = (E(1,:) > 0) & (E(2,:) > 0);
plot3(x(1,i),x(2,i),x(3,i),'or'), hold on
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
title('Flow','Fontsize',16), axis square

return



%% subroutines
%==========================================================================

function F = spm_ness_F(P,M,U,OPT)
% equations of motion for flow modelling
% P - model parameters
% M - gmodel structure
% U - design or basis functions
%__________________________________________________________________________
if nargin > 3

    F = spm_NESS_gen_lap(P,M,U,OPT);

else

    % autonomous flow
    %--------------------------------------------------------------------------
    F = spm_NESS_gen_lap(P,M,U);

end

function f = spm_ness_f(x,u,P,M,OPT)
% equations of motion for DEM (filtering)
% x - latent state
% u - exogenous input
% P - model parameters
% M - model structure
%__________________________________________________________________________
if nargin > 4
    f = spm_fx_NESS(x,u,P,M,OPT);
    return
elseif nargin > 3
    f = spm_fx_NESS(x,u,P,M);
else
    f = spm_fx_NESS(x,u,P);
end


return

function X = spm_grad(X)
% gradients based on diff operator (i.e. the past)
% X  - matrix of variables
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
n  = size(U,1);
X  = spm_dctmtx(n,d);
B  = X(T,:)\U(T,:);
U  = U - X*B;

% get eigenvectors
%--------------------------------------------------------------------------
u  = U(T,:);
v  = spm_svd(spm_en(u)');
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


%% Iterated simulations
%==========================================================================
spm_clear, clear
n     = 10;
tab   = zeros(5,4,n);
DEM   = cell(1,n);
F     = cell(1,n);
for i = 1:n

    % durations
    %----------------------------------------------------------------------
    SIM.N  = (i - 1)*52;        % end point (weeks)
    SIM.D  = 256;               % depth of training data (weeks)
    SIM.nT = 52;                % duration of simulation (in weeks)
    SIM.dT = 2;                 % time between rebalancing (in weeks)

    % specify number of indicator states and assets
    %----------------------------------------------------------------------
    SIM.n  = 3;                 % number of indicator states
    SIM.m  = 6;                 % number of assets
    SIM.d  = 8;                 % order of detrending
 
    [Tab,f,Dem] = DEM_FIN(SIM);
    tab(:,:,i) = table2array(Tab);
    DEM{i}     = Dem; 
    F{i}       = f;
end

% bar chart results
%--------------------------------------------------------------------------
spm_figure('GetWin','Annual performance'); clf
for i = 1:4
    subplot(3,2,i)
    bar(squeeze(tab(:,i,:))')
    title(Tab.Properties.VariableNames{i})
    legend(Tab.Properties.RowNames)
    axis square
end

L = full(spm_cat(F'));

subplot(3,3,7)
bar(L(:,1)), xlabel('year'),ylabel('nats')
title('ELBO')
axis square

subplot(3,3,8)
bar(-L(:,2)), xlabel('year'),ylabel('nats')
title('Coupling')
axis square

subplot(3,3,9)
bar(-L(:,3)), xlabel('year'),ylabel('nats')
title('Chaos')
axis square

avg  = mean(tab,3);
VariableNames{1} = 'Annual RoR (%)';
VariableNames{2} = 'Volatility (%)';
VariableNames{3} = 'Sharpe ratio';
VariableNames{4} = 'Drawdown (%)';
RowNames = {'hold','max-ret','risk-ret','max-pro','risk-pro'};

Avg = array2table(avg);
Avg.Properties.VariableNames = VariableNames;
Avg.Properties.RowNames      = RowNames


save DEMFIN


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
