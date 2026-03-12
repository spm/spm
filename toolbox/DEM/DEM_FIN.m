function [Tab,F,DEM] = DEM_FIN(SIM,OPT)
% FORMAT [Tab,F,DEM] = DEM_FIN(SIM)
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
    D    = 8*32;                % depth of training data (weeks)
    nT   = 52;                  % duration of simulation (in weeks)
    dT   = 4;                   % time between rebalancing (in weeks)

    % specify number of indicator states and assets
    %----------------------------------------------------------------------
    n    = 3;                   % number of indicator states
    m    = 6;                   % number of assets
    d    = 2;                   % order of detrending

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

% exogenous inputs 
%--------------------------------------------------------------------------
U    = spm_dctmtx(D,d,1:(t + 32));
U    = U*sqrt(D);

% Exchange-Traded Funds
%--------------------------------------------------------------------------
V    = spm_log(V(:,1:m));       % log value
L    = spm_grad(V);             % rate of log return per day
V    = K*V;                     % log return per day
L    = K*L*dt;                  % rate of log return per week

% detrend and eigenreduce log indicators
%--------------------------------------------------------------------------
O    = spm_log(O);              % log indicators
O    = K*O;                     % decimate
t    = (T - D):T;               % past period
I    = spm_eigenreduce([O,L],t);
I    = I(:,1:n);

% plot data
%--------------------------------------------------------------------------
spm_figure('GetWin','Data'); clf;

% response and explanatory variables
%--------------------------------------------------------------------------
f    = (T - D):(T + nT);        % total period
subplot(3,1,1), hold off
plot(date(f),I(f,1:n)), title('Log indicator variables','FontSize',14)
datetick('x','mmm-yy'), xlabel('time (mmm-yy'), spm_axis tight

subplot(3,1,2), hold off
plot(f,V(f,:)), title('Log value of ETFs','FontSize',14)
xlabel('time (days)'), spm_axis tight

subplot(3,1,3), hold off
plot(f,L(f,:)/sqrtm(cov(L(t,:)))), hold on
plot(f,2*U(f,:),'--')
plot([T,T],    get(gca,'YLim'),'-.r')

title('Log return (per week)','FontSize',14)
xlabel('time (weeks)'), spm_axis tight


%% system identification
%==========================================================================

% detrend and normalise
%--------------------------------------------------------------------------
gx  = @(x,u,P) P.scale'*x + P.trend'*u;
gy  = @(y,u,P) (y - u*P.trend)/P.scale;

u   = [U(t,:)];                 % exogenous inputs 
Y   = [I(t,:) L(t,:)];          % response variables

P.trend = u\Y;
P.scale = sqrtm(cov(Y - u*P.trend));

% latent states
%--------------------------------------------------------------------------
X   = gy(Y,u,P);

% prior over precision of states and random fluctuations
%--------------------------------------------------------------------------
C   = cov(X);
W   = cov(gradient(X')');
W   = diag(W);
W   = W/4;
W   = inv(diag(W));

% get priors over model parameters (polynomial coeficients)
%--------------------------------------------------------------------------
K     = 2;                      % order of polynomial expansion plus one
P.Qp  = 1/32;                   % variance of parameters (solenoidal)
P.Sp  = 0;                      % variance of parameters (surprisal)
P.Rp  = 1;                      % variance of parameters (expectation)

% constraints on solenoidal flow
%--------------------------------------------------------------------------
J      = cell(3,3);
J{1,1} = ones(n,n);
J{2,2} = ones(m,m);
J{2,1} = ones(m,n);
J{1,2} = ones(n,m);
J      = full(spm_cat(J));

% cconstraints on surprisal flow
%--------------------------------------------------------------------------
Q      = cell(3,3);
Q{1,1} = zeros(n,n);
Q{2,2} = ones(m,m);
Q{2,1} = ones(m,n);
Q{1,2} = zeros(n,m);
Q      = full(spm_cat(Q));

% preclude states enslaving themselves for long term forecasting
%--------------------------------------------------------------------------
if ~nargin
    Q = Q - diag(diag(Q));
end

% get priors 
%--------------------------------------------------------------------------
[pE,pC] = spm_NESS_priors(size(J,1),K,P,W,J,Q,C);

% mean
%--------------------------------------------------------------------------
pE.Rp(1,:) = 0;
pC.Rp(1,:) = 0;

% precision
%--------------------------------------------------------------------------
pC.Sp(1,:,:) = 0;

% state-dependent fluctuations
%--------------------------------------------------------------------------
pC.Gp = Q/32;

% scaling
%--------------------------------------------------------------------------
pE.trend = P.trend;
pC.trend = P.trend*0;
pE.scale = P.scale;
pC.scale = P.scale*0;

%% model inversion with Variational Laplace (fitting flow)
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
M.Nmax = 32;                   % maximum number of iterations
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
f      = gradient(X')';        % target flow

% posterior over parameters
%--------------------------------------------------------------------------
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,B,f);


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
% [Ep,Cp,Eh,F] = spm_nlsi_GN(M,B,f);
%--------------------------------------------------------------------------
% or iterate if necessary
%--------------------------------------------------------------------------
% M.hE         = Eh;
% M.pE.W       = diag(exp(Eh));
% [Ep,Cp,Eh,F] = spm_nlsi_GN(M,B,f);
%--------------------------------------------------------------------------


% posterior precision of random fluctuations
%--------------------------------------------------------------------------
F0     = F;                    % ELBO     
W      = diag(exp(Eh));        % posterior volatility
Ep.W   = W;                    % posterior volatility
M.W    = W;                    % update for flow model

% Bayesian belief updating
%--------------------------------------------------------------------------
M.pE   = Ep;                   % and parameters of flow
M.pC   = Cp;                   % and covariance of parameters
M.hE   = Eh;                   % and parameters of flow
M.hC   = 1/64;                 % and covariance of parameters


%% evaluate Jacobian at initial state
%==========================================================================
spm_figure('GetWin','Jacobian'); clf;

% numerical Jacobian
%--------------------------------------------------------------------------
spm_J = @(x,u,P) full(spm_diff(@spm_fx_NESS,x,u,P,[],1));
J     = spm_J(x0,u0,Ep);

fprintf('\nJacobian at initial conditions\n')
disp(J)
subplot(3,2,1)
imagesc(J), title('Jacobian (causal coupling)','FontSize',14)
xlabel('latent states'),ylabel('latent states')
axis square

% Lyapunov exponents
%--------------------------------------------------------------------------
[e,v] = eig(J);
v     = diag(v);
[~,i] = sort(real(v),'descend');
v     = v(i);
e     = e(:,i);

% plot
%--------------------------------------------------------------------------
str   = {};
for i = 1:n, str = [str, 'Ind']; end
str   = [str, EFT];

subplot(3,2,2), plot(v,'.k','MarkerSize',16)
hold on, plot([0,0],get(gca,'YLim'),':k'), hold off
hold on, plot(get(gca,'XLim'),[0 0],':k'), hold off
title('Lyapunov exponents','FontSize',14)
xlabel('real part'),ylabel('imaginary part')
axis square
subplot(3,2,3), bar(-1./real(v),'m')
hold on, plot([1,numel(v)],[4 4],'--k'), hold off
title('Time constants','FontSize',14)
xlabel('eigenmode'),ylabel('weeks')
axis square
subplot(3,2,4), bar(imag(v)/(2*pi))
title('Frequency','FontSize',14)
xlabel('eigenmode'),ylabel('cycles per week')
axis square
subplot(3,1,3), bar(abs(e'))
title('Eigenmodes','FontSize',14)
xlabel('eigenmode'),ylabel('abolute value'), legend(str)

% Lyapunov exponent and Hausdorff dimension (Kaplan-Yorke conjecture)
%--------------------------------------------------------------------------
ni    = 64;
nx    = n + m;
E     = zeros(nx,ni);
F     = zeros(nx,ni);
Z     = zeros(nx,ni);
J     = zeros(nx,nx,ni);
for i = 1:ni
    F(:,i)   = spm_fx_NESS(X(i,:),U(i,:),Ep);
    J(:,:,i) = spm_J(X(i,:),U(i,:),Ep);
    Z(:,i)   = eig(J(:,:,i));
    E(:,i)   = sort(Z(:,i),'descend','ComparisonMethod','real');
end

% Lyapunov exponent and Hausdorff dimension (Kaplan-Yorke conjecture)
%--------------------------------------------------------------------------
LE    = mean(real(E),2);
j     = sum(LE >= 0);
HD    = j + sum(LE(1:j))/abs(LE(j + 1));
V0    = min(LE);
fprintf('Hausdorff Dimension: %.2f\n',HD)


% serial correlations
%==========================================================================
if nargin

    % smoothness
    %----------------------------------------------------------------------
    s = 1/2;

else
    
    % evaluate from residuals of flow
    %----------------------------------------------------------------------
    spm_figure('GetWin','Serial correlations'); clf;

    % residuals of flow
    %----------------------------------------------------------------------
    r = f - spm_ness_F(Ep,M,B);
    r = r/diag(std(r));
    subplot(2,2,1), hold off
    for i = 1:n
        plot(xcorr(r(:,i),32,'normalized')), hold on
    end
    xlabel('lag'), ylabel('correlation'), axis square, title('indicators')
    subplot(2,2,2), hold off
    for i = (1:m) + n
        plot(xcorr(r(:,i),32,'normalized')), hold on
    end
    xlabel('lag'), ylabel('correlation'), axis square, title('EFTs')

    % serial correlations and smoothness
    %----------------------------------------------------------------------
    s     = r(:,end);
    for i = 1:2
        s = [s gradient(s(:,end)')'];
    end
    subplot(2,2,3), rho = cov(s); imagesc(rho)
    xlabel('order of motion'), axis square, title('Covariance: empirical')
    s     = sqrt((1/2)*1/rho(2,2));
    [~,c] = spm_DEM_R(3,s);
    subplot(2,2,4), imagesc(c)
    xlabel('order of motion'), axis square, title('Covariance: analytic')

    fprintf('\nsmoothness (s) = %.2f\n',s);

end


%% Bayesian model comparison
%==========================================================================

% get reduced priors: no enlaving
%--------------------------------------------------------------------------
J        = cell(3,3);
J{1,1}   = ones(n,n);
J{2,2}   = ones(m,m);
J        = spm_cat(J);

% get priors 
%--------------------------------------------------------------------------
[rE,rC]  = spm_NESS_priors(size(J,1),K,P,W,J);
rE.scale = pE.scale;
rC.scale = pC.scale;
rE.trend = pE.trend;
rC.trend = pC.trend;

F1       = spm_log_evidence_reduce(Ep,Cp,pE,pC,rE,rC);
fprintf('\nLog-evidence for enslaving: %.2f\n',F1)

% get reduced priors: predicability of indicator variables
%--------------------------------------------------------------------------
J        = cell(3,3);
J{1,1}   = eye(n,n);   %%% remove coupling
J{2,2}   = eye(m,m);   %%% remove coupling
J        = spm_cat(J);

% get priors 
%--------------------------------------------------------------------------
[rE,rC]  = spm_NESS_priors(size(J,1),K,P,W,J);
rE.trend = pE.trend;
rC.trend = pC.trend;
rE.scale = pE.scale;
rC.scale = pC.scale;

F2       = spm_log_evidence_reduce(Ep,Cp,pE,pC,rE,rC);
fprintf('\nLog-evidence for solenoidal coupling: %.2f\n',F2)



%% model inversion with generlized filtering
%==========================================================================
DEM.EFT = M.EFT;
DEM.G   = M;                   % generative model of flow


% serial correlations and orders of motion
%--------------------------------------------------------------------------
DEM.M(1).E.s      = s;         % smoothness of fluctuations
DEM.M(1).E.n      = 1;         % order of generalised motion (states)
DEM.M(1).E.d      = 1;         % order of generalised motion (data)
DEM.M(1).E.nD     = 2;         % number of integrations per dt
DEM.M(1).E.nE     = 4;         % number of iterations
DEM.M(1).E.nN     = 4;         % number of iterations
DEM.M(1).E.linear = 4;         % differentiation scheme
DEM.M(1).E.v      = V0;        % lyapunov exponent

% model
%--------------------------------------------------------------------------
DEM.M(1).f  = @spm_fx_NESS;     % flow
DEM.M(1).g  = gx;              % observer function
DEM.M(1).gy = gy;              % observer function
DEM.M(1).x  = x0;              % initial state
DEM.M(1).pE = Ep;              % posterior esimates from VL
DEM.M(1).pC = Cp;              % posterior esimates from VL
DEM.M(1).V  = exp(16);
DEM.M(1).W  = W;

DEM.M(2).x  = [];
DEM.M(2).v  = u0;              % initial cause
DEM.M(2).V  = exp(16);         % prior precicion


% Invert using generalised filtering: in generalised coordinates of motion 
%--------------------------------------------------------------------------
DEM.Y = Y';
DEM.U = u';

% DEM   = spm_LAP(DEM);
% 
% % illustrate results
% %------------------------------------------------------------------------
% spm_DEM_qU(DEM.qU)
% subplot(2,2,1), hold on, plot(DEM.Y','.k'), hold off


% return modelling results
%--------------------------------------------------------------------------
F   = [F0,F1,F2,V0];

if nargin > 1

    % dummy forecasting table
    %----------------------------------------------------------------------
    Tab = array2table(zeros(5,4));
    return
    
end


%% forecasting - analytic
%==========================================================================
spm_figure('GetWin','Fokker Planck'); clf

% future causes
%--------------------------------------------------------------------------
Nf  = 32;                      % forecasting period
t   = (1:Nf) + T;              % future

DEM.U      = U(t,:)';          % future causes
DEM.X      = X';               % past consquences
DEM.Y      = Y';               % past consquences
DEM.M(1).x = X(end,:)';        % initial state
DEM.M(2).v = u(end,:)';        % initial cause

[Ez,Cz,Vz] = spm_NESS_forecast(DEM);

% overlay actual outcomes if specified
%--------------------------------------------------------------------------
xLim  = [(D - Nf),(D + Nf)];
s     = (1:Nf) + D;
for i = 1:n
    subplot(6,2,i), hold on, set(gca,'ColorOrderIndex',i),
    plot(s + 1,I(t,i),'.')
    set(gca,'XLim',xLim)
end

% enslaved outcomes
%--------------------------------------------------------------------------
for i = (1:m) + n
    subplot(6,2,i), hold on, set(gca,'ColorOrderIndex',i),
    plot(s + 1,L(t,i - n),'.')
    set(gca,'XLim',xLim)
end

%% forecasting - numerical
%==========================================================================
if ~nargin

    % predictive posterior realizations
    %----------------------------------------------------------------------
    spm_figure('GetWin','Forecasting'); clf
    [Ez,Cz,Vz,Py] = spm_NESS_forecasting(DEM);

    % overlay actual outcomes if specified
    %----------------------------------------------------------------------
    subplot(2,2,1), hold on, set(gca,'ColorOrderIndex',1),
    for i = 1:n
        plot(s,I(t,i),'.')
    end

    % overlay actual outcomes if specified
    %----------------------------------------------------------------------
    xLim  = [(D - 4*Nf),(D + Nf)];
    for i = 1:n
        subplot(8,3,12 + i), hold on, set(gca,'ColorOrderIndex',i),
        plot(s,I(t,i),'.')
        set(gca,'XLim',xLim)
    end

    % enslaved outcomes
    %----------------------------------------------------------------------
    for i = (1:m) + n
        subplot(8,3,12 + i), hold on, set(gca,'ColorOrderIndex',i),
        plot(s,L(t,i - n),'.')
        set(gca,'XLim',xLim)
    end

    % Projections
    %----------------------------------------------------------------------
    spm_figure('GetWin','Responses'); clf

    % extract response variables
    %----------------------------------------------------------------------
    z     = (1:nT) + T;            % future time points with data
    t     = (1:Nf) + T - 1;        % future time points predicted
    i     = (1:m) + n;             % indices of response variables
    EL    = Ez(i,:);               % posterior predictive expectations
    VL    = Vz(i,:);               % posterior predictive variances

    % plot credible intervals and data
    %----------------------------------------------------------------------
    xLim  = [T - 2*Nf,T + Nf];
    for i = 1:m

        % predictive posterior
        %------------------------------------------------------------------
        subplot(6,2,i), hold off
        set(gca,'ColorOrderIndex',i + n)
        spm_plot_ci(EL(i,:),VL(i,:),t), hold on

        % overlay past and future outcomes
        %------------------------------------------------------------------
        set(gca,'ColorOrderIndex',i + n)
        plot(1:T,L(1:T,i),'LineWidth',2), hold on
        set(gca,'ColorOrderIndex',i + n)
        plot(z,L(z,i),'.'), plot(xLim,[0,0],':'),hold off
        title(['RoR: ' EFT{i}],'FontSize',14)
        set(gca,'XLim',xLim)
    end
end


%% simulate portfolio management under this generative model
%==========================================================================
if nargin, DEM.nograph = 1; end

% simulate portfolio management
%--------------------------------------------------------------------------
Tab = spm_portfolio(L,I,U,DEM,T,dT)

if nargin, return, end




% free energy and Market cycles
%==========================================================================
spm_figure('GetWin','Financial free energy'); clf

% surprisal in the past
%--------------------------------------------------------------------------
F     = [];
d     = 52*3;
for t = (T - d):(T + Nf - 1)
    u    = U(t,:);
    x    = gy([I(t,:) L(t,:)],u,Ep);
    F(t) = spm_fx_NESS(x,u,Ep,[],'S');
end

% surprisal forecasts
%--------------------------------------------------------------------------
Np    = size(Py,3);
f     = zeros(T + Nf - 1,Np);
for i = 1:Np
    for t = 1:Nf
        u              = U(T + t,:);
        x              = gy(Py(:,t,i)',u,Ep);
        f(T + t - 1,i) = spm_fx_NESS(x,u,Ep,[],'S');
    end
end

% free energy fluctuations
%--------------------------------------------------------------------------
subplot(2,1,1), hold off
t = M.date;
i = (T - 0):(T + Nf - 1);
plot(t(i),f(i,:),'r:','LineWidth',1), hold on
i = (T - d):(T + Nf - 1);
plot(t(i),F(i),'b','LineWidth',1), set(gca,'XLim',[t(i(1)),t(i(end))])
plot(get(gca,'XLim'),[5,5],'-.k')
datetick('x','mmm-yy'), xlabel('time'), ylabel('surprisal (nats)');
title('Financial free energy','FontSize',14)

% Phase portraits: interpolate free energy fluctuations
%--------------------------------------------------------------------------
S     = interp(F,8);
for i = 1:Np
    Si(:,i) = interp(f(:,i),8);
end

subplot(2,1,2), hold off
dSdt  = gradient(S);
dFdt  = gradient(Si')';
plot(dFdt,Si,'r:','LineWidth',1), hold on
plot(dSdt,S,'b','LineWidth',1)
plot([0,0],get(gca,'YLim'),':k')
axis square, title('Phase portrait','FontSize',14)
xlabel('time derivative'), ylabel('surprisal (nats)');


%% Is the Market chaotic?
%==========================================================================
spm_figure('GetWin','Flow'); clf

% state space for evaluation
%--------------------------------------------------------------------------
ni    = 64;
nx    = n + m;
E     = zeros(nx,ni);
F     = zeros(nx,ni);
Z     = zeros(nx,ni);
J     = zeros(nx,nx,ni);
for i = 1:ni
    F(:,i)   = spm_fx_NESS(X(i,:),U(i,:),Ep);
    J(:,:,i) = spm_J(X(i,:),U(i,:),Ep);
    Z(:,i)   = eig(J(:,:,i));
    E(:,i)   = sort(real(Z(:,i)),'descend');
end

subplot(2,2,1)
imagesc(mean(J,3)), title('Jacobian (causal coupling)','FontSize',14)
xlabel('latent states'),ylabel('latent states'), axis square

subplot(2,2,2), plot(Z,'b.','MarkerSize',8)
hold on, plot([0,0],get(gca,'YLim'),':k'), hold off
hold on, plot(get(gca,'XLim'),[0 0],':k'), hold off
title('Lyapunov exponents','FontSize',14)
xlabel('real part'),ylabel('imaginary part')
axis square


% generate path of least action from initial conditions
%-------------------------------------------------------------------------
G      = DEM.M;         % generative model
G(1).E.nD = 4;          % number of integrations per dt

ni     = 128;
t      = (T - ni):T;
G(1).x = x0;
G(1).V = exp(16);
G(1).W = exp(16);
G      = spm_DEM_generate(G,U(t,:)');

x = G.pU.x{1};
X = x';
E = zeros(nx,ni);
F = zeros(nx,ni);
Z = zeros(nx,ni);
J = zeros(nx,nx,ni);
for i = 1:ni
    F(:,i)   = spm_fx_NESS(X(i,:),U(i,:),Ep);
    J(:,:,i) = spm_J(X(i,:),U(i,:),Ep);
    Z(:,i)   = eig(J(:,:,i));
    E(:,i)   = sort(real(Z(:,i)),'descend');
end

subplot(2,1,2), hold off
i     = 1:ni;
quiver3(x(1,i),x(2,i),x(3,i),F(1,i),F(2,i),F(3,i),'b'), hold on
plot3(x(1,i),x(2,i),x(3,i),'.b','MarkerSize',8)
plot3(x(1,i),x(2,i),x(3,i),':b')
i     = (E(1,:) > 0) & (E(2,:) > 0);
plot3(x(1,i),x(2,i),x(3,i),'.r','MarkerSize',16), hold on
quiver3(x(7,i),x(8,i),x(9,i),F(7,i),F(8,i),F(9,i),'m'), hold on
plot3(x(7,i),x(8,i),x(9,i),'.m','MarkerSize',8)
plot3(x(7,i),x(8,i),x(9,i),':m')
i     = (E(1,:) > 0) & (E(2,:) > 0);
plot3(x(7,i),x(8,i),x(9,i),'.r','MarkerSize',16), hold on
plot3(0,0,0,'.g','MarkerSize',32), hold on
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
title('Paths of least action','Fontsize',16), axis square


return



%% subroutines
%==========================================================================

function f = spm_ness_F(P,M,U,OPT)
% equations of motion for flow modelling
% P - model parameters
% M - gmodel structure
% U - design or basis functions
%__________________________________________________________________________
if nargin > 3
    f = spm_NESS_gen_lap(P,M,U,OPT);
else
    f = spm_NESS_gen_lap(P,M,U);
end

% function f = spm_fx_NESS(x,u,P,M,OPT)
% equations of motion for DEM (filtering)
% x - latent state
% u - exogenous input
% P - model parameters
% M - model structure
%__________________________________________________________________________


return

function X = spm_grad(X)
% gradients based on diff operator (i.e. the past)
% X  - matrix of variables
%__________________________________________________________________________

X  = diff(X);
X  = [X(1,:); X];

return


function I = spm_eigenreduce(U,t,L)
% eigenreduction of indicator variables
% FORMAT I = spm_eigenreduce(U,t,L)
%--------------------------------------------------------------------------
% U  - indicator variables
% t  - time points
% L  - response variables
%
% I  - normalised canonical variates
% W  - canonical vectors (weights)
%__________________________________________________________________________

% detrend
%--------------------------------------------------------------------------
X  = spm_dctmtx(size(U,1),2);
U  = U - X*(X(t,:)\U(t,:));

if nargin < 3

    % singular value decomposition
    %----------------------------------------------------------------------
    u   = spm_en(U(t,:));
    v   = spm_svd(u');

else

    % canonical variate analysis
    %----------------------------------------------------------------------
    CVA = spm_cva(gradient(L(t,:)')',U(t,:));
    v   = CVA.W;

end

% normalise
%--------------------------------------------------------------------------
I  = U*v;
I  = I/diag(std(I));

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
J{2,2} = ones(m,m);
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
cell2mat(spm_NESS_gen_lap(Ep,M,x,'DS'))

full(spm_diff(@spm_fx_NESS,x,[],Ep,[],'S',1)')
full(spm_diff(@spm_NESS_gen_lap,Ep,M,x,'S',3)')

% effect of parameters on Jacobian
%--------------------------------------------------------------------------
dfdP  = spm_diff(@spm_fx_NESS,x,[],Ep,[],3);
dJdP  = spm_diff(@spm_fx_NESS,x,[],Ep,[],[1 3]);

subplot(2,1,1)
imagesc(dfdP)






