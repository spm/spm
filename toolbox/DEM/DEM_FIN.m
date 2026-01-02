function DCM = DEM_FIN
% FORMAT DCM = DEM_FIN
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
rng(1)


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

i    = [2 6 3 5 4 1];
V    = V(:,i);
EFT  = EFT(i);


% remove null indicator variables
%--------------------------------------------------------------------------
d    = find(min(U) > 0);
U    = U(:,d);
name = name(d);

% decimate and smooth
%--------------------------------------------------------------------------
dt   = 5;                       % sampling interval (e.g., a week)
t    = numel(date);             % number of samples
d    = 1:dt:t;                  % time series for eigenreduction
nT   = 64;                      % duration of simulation (in dt)
dT   = 4;                       % time between rebalancing (in dt)
T    = numel(d) - nT;           % start of simulation
D    = 512;                     % depth of training data
date = date(d);

% Exchange-Traded Funds
%--------------------------------------------------------------------------
V    = spm_log(V);              % log value
L    = gradient(V')';           % rate of log return per day
L    = spm_conv(L,dt,0);        % smooth before decimating
L    = L(d,:)*dt;               % decimated RoR per dt

% detrend in log space to condition data under NESS assumption
%--------------------------------------------------------------------------
U    = spm_log(U);              % log indicators
U    = spm_conv(U,16,0);      % smoothing
U    = U(d,:);                  % decimated indicators
I    = spm_eigenreduce(U,1:T,L);% eigenvariates


% specify
%--------------------------------------------------------------------------
n    = 4;                     % number of indicator states
m    = 6;                     % number of assets

% plot data
%--------------------------------------------------------------------------
spm_figure('GetWin','Data'); clf;

% response and explanatory variables
%--------------------------------------------------------------------------
t    = (T - D):T;
I    = [I(:,1:n)];            % indicator variables
L    = [L(:,1:m)];            % rate of log return
Y    = [I(t,:),L(t,:)];       % training data

subplot(2,1,1)
plot(date,I), title('Log indicator variables','FontSize',14)
datetick('x','mmm-yy'), xlabel('time')

subplot(2,1,2)
plot(L(:,1:m)), title('Log return (%)','FontSize',14)
xlabel('time')


% system identification
%==========================================================================

% scaling
%--------------------------------------------------------------------------
scale = std(Y);
X     = rdivide(Y,scale);

% prior over precision of states and random fluctuations
%--------------------------------------------------------------------------
W     = var(gradient(X')');
W     = W/(8^2);
W     = diag(1./W);

% model inversion with Variational Laplace (fitting flow)
%==========================================================================

% get model parameters (polynomial coeficients)
%--------------------------------------------------------------------------
K     = 2;                    % order of polynomial expansion plus one
P     = 1/32;                 % precision of parameters (prior)

% constraints on model parameters (polynomial coefficients)
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
Q{2,1} = ones(m,n);
Q      = spm_cat(Q);

% get priors 
%--------------------------------------------------------------------------
[pE,pC] = spm_NESS_priors(size(J,1),K,P,W,J,Q);

% scaling
%--------------------------------------------------------------------------
pE.scale = scale(:);
pC.scale = pE.scale*0;

%% Variational Laplace
%--------------------------------------------------------------------------

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

% update initial expectations (parameters)
%----------------------------------------------------------------------
M.P    = Ep;

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
disp(full(spm_diff(@spm_fx_NESS,x0,[],Ep,[],1)))

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
fprintf('\nLog-evidence for stochastic choas: %.2f\n',F)



% model inversion with generlized filtering
%==========================================================================

% serial correlations and orders of motion
%--------------------------------------------------------------------------
DEM.M(1).E.s      = 0;         % smoothness of fluctuations
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
% %--------------------------------------------------------------------------
% spm_DEM_qU(DEM.qU)
% subplot(2,2,1), hold on, plot(DEM.Y','.k'), hold off



%% forecasting - analytic
%==========================================================================
spm_figure('GetWin','Fokker Planck'); clf

Nf         = max(2*dT,8);      % forecasting period
DEM.T      = T;
DEM.Y      = Y';               % legacy (training) data
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
t     = (1:Nf) + T;
for i = 1:n
    subplot(6,2,i), hold on, set(gca,'ColorOrderIndex',i),
    plot(t,I(t,i),'.')
end

% enslaved outcomes
%--------------------------------------------------------------------------
for i = (1:m) + n
    subplot(6,2,i), hold on, set(gca,'ColorOrderIndex',i),
    plot(t,L(t,i - n),'.')
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

% enslaving outcomes
%--------------------------------------------------------------------------
subplot(4,1,3), hold on, set(gca,'ColorOrderIndex',1),
for i = 1:n
    plot(t,I(t,i),'.')
end

% enslaved outcomes
%--------------------------------------------------------------------------
for i = 1:m
    subplot(4,1,4), hold on
    set(gca,'ColorOrderIndex',i),
    plot(t,L(t,i),'.')
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

% free energy and Market cycles
%==========================================================================

% NESS density and expected flow
%--------------------------------------------------------------------------
% [F,S] = spm_NESS_gen_lap(Ep,M);



% simulate portfolio management under this genertive model
%==========================================================================
spm_figure('GetWin','Portfolio forecasts'); clf

M.Nmax = 32; % maximum number of iterations
tab    = spm_portfolio(L,I,M,DEM,T,dT);

return


function I = spm_eigenreduce(U,T,L)
% eigenreduction of indicator variables
% U  - indicator variables
% T  - based on first T samples
% L  - response variables
%__________________________________________________________________________
n    = size(U,1);
X    = spm_dctmtx(n,4);
B    = X(T,:)\U(T,:);
u    = U(T,:) - X(T,:)*B;
U    = U - X*B;

% get eigenvectors
%--------------------------------------------------------------------------
v    = spm_svd(spm_en(u)');
I    = U*v;
I    = spm_en(I)*sqrt(n);

return

% Augment explanatory and response variables
%--------------------------------------------------------------------------
n    = size(U,1);
X    = spm_en([U, spm_dctmtx(n,4)]);

% Canonical variates analysis
%--------------------------------------------------------------------------
s    = numel(T)/16;
Y    = spm_conv(L,s,0);
Y    = spm_detrend(Y(T,:));
CVA  = spm_cva(Y,X(T,:),[],[],8);
I    = X*CVA.W;
I    = spm_en(I)*sqrt(n);


return

function tab = spm_portfolio(L,I,M,DEM,Ti,dT)
% FORMAT spm_portfolio(L,DEM)
% simulate portfolio managment
%--------------------------------------------------------------------------
% L   - rate of log return per M.dt
% I   - indicator variables
% M   - flow model
% DEM - generative model
% Ti  - time of initial investment
% dT  - time between rebalancing
%__________________________________________________________________________

% policies
%--------------------------------------------------------------------------
policies = {'Buy and hold', ...
    'Retrospective: Expected RoR'...
    'Retrospective: Risk sensitive'...
    'Prospective: Expected RoR'...
    'Prospective: Risk sensitive'...
    };

% get sizes and times
%--------------------------------------------------------------------------
dt    = M.dt;                 % period (days)
Ann   = 100*365/dt;           % scaling for annualised % RoR
Tf    = size(L,1);            % time of final investment
m     = size(L,2);            % number of (aggregated) assets

% prior preferences
%--------------------------------------------------------------------------
try  RoR    = DEM.RoR;    catch, RoR    = 100; end
try  Sharpe = DEM.Sharpe; catch, Sharpe = 1;  end

m_p   = RoR/Ann;              % prior RoR per dt
s_p   = m_p/Sharpe;           % prior volitility (s.d.)
c_p   = s_p^2;                % prior volitility (variance)

% set up policies: moving funds from one assets to another nP times
%--------------------------------------------------------------------------
dP    = .5;                   % proportional transaction
cP    = .001;                 % proportional cost of transaction
nP    = 3;                    % number of transactions
p     = spm_combinations([m,m]);
np    = size(p,1);
Pk    = cell(1,np);
for k = 1:np
    i    = p(k,1);
    j    = p(k,2);
    Pij  = eye(m,m);

    % transaction (from asset j to i)
    %----------------------------------------------------------------------
    Pij(i,j) = Pij(i,j) + dP;
    Pij(j,j) = Pij(j,j) - dP;

    % cost of transaction
    %----------------------------------------------------------------------
    if i ~= j
        Pij(i,j) = Pij(i,j) - cP;
    end

    % k-th transaction
    %----------------------------------------------------------------------
    Pk{k} = Pij;
end

% compound policies; i.e., compositions of unitary transactions
%--------------------------------------------------------------------------
p     = spm_combinations(ones(1,nP)*np);
np    = size(p,1);
P     = cell(1,np);
for i = 1:size(p,1)
    Pij   = 1;
    for j = 1:size(p,2)
        Pij = Pij*Pk{p(i,j)};
    end
    P{i}  = Pij;
end

% simulate allocations and ensuing resturns
%==========================================================================
Np    = 5;                     % number of policies to simulate
W     = ones(m,Np)/m;          % inital (flat) weights
R     = zeros(Tf,Np);          % cumulative RoR
D     = zeros(Tf,Np);          % RoR per dt
l     = (1:m) + size(I,2);     % indices of RoR
for i = Ti:Tf

   
   % cumulative return on investment per dt
   %-----------------------------------------------------------------------
   D(i,:) = L(i,1:m)*W;
   R(i,:) = R(i - 1,:) + D(i,:);

   % portfolio allocations
   %-----------------------------------------------------------------------
   Wret(:,i)  = W(:,3);
   Wpro(:,i)  = W(:,5);

    % transactions
    %-----------------------------------------------------------------------
    if ~rem(i,dT)

        % buy and hold 
        %------------------------------------------------------------------
        W(:,1)  = P{1}*W(:,1);

        % retrospective
        %==================================================================
        ti    = i - (1:dT) + 1;
        EL    = L(i,1:m);
        CL    = cov(L(ti,1:m));
        F     = zeros(1,np);
        G     = zeros(1,np);
        for k = 1:np

            % expected utility
            %--------------------------------------------------------------
            W_k  = P{k}*W(:,2);
            m_q  = EL*W_k;
            F(k) = m_q;

            % risk sensitive
            %--------------------------------------------------------------
            W_k  = P{k}*W(:,3);
            m_q  = EL*W_k;
            c_q  = W_k'*CL*W_k;
            G(k) = spm_kl_normal(m_q,c_q,m_p,c_p);

        end

        % expected utility
        %------------------------------------------------------------------
        [~,j]  = max(F);
        W(:,2) = P{j}*W(:,2);

        % risk sensitive
        %------------------------------------------------------------------
        [~,j]  = min(G);
        W(:,3) = P{j}*W(:,3);


        % prospective
        %==================================================================
        if false

            % upper bound on performance, if the future were known
            %--------------------------------------------------------------
            for t = 1:dT
                try
                    Ey(t,:) = L(i + t,1:m);
                    Cy{t}   = cov(L((1:dT) + i,1:m));
                catch
                    Ey(t,:) = L(end,1:m);
                    Cy{t}   = cov(L(end - (1:dT) + 1,1:m));
                end
            end
        else

            % posterior predictive density over RoR
            %--------------------------------------------------------------
            j          = 1:i;
            [Ez,Cz,M]  = spm_forecast_update(L(j,:),I(j,:),M,DEM,dT + 1);
            for t = 1:dT
                Ey(t,:) = Ez(l,t + 1)';
                Cy{t}   = Cz{t + 1}(l,l);
            end

        end

        % for each policy
        %------------------------------------------------------------------
        F     = zeros(1,np);
        G     = zeros(1,np);
        for k = 1:np

            % expected utility
            %--------------------------------------------------------------
            W_k   = P{k}*W(:,4);
            for t = 1:dT

                % predictive posterior over outcomes
                %----------------------------------------------------------
                m_q  = Ey(t,:)*W_k;

                % path integral of expected free energy (risk)
                %----------------------------------------------------------
                F(k) = F(k) + m_q;

            end

            % risk sensitive
            %--------------------------------------------------------------
            W_k   = P{k}*W(:,5);
            for t = 1:dT

                % predictive posterio over outcomes
                %----------------------------------------------------------
                m_q  = Ey(t,:)*W_k;
                c_q  = W_k'*Cy{t}*W_k;

                % path integral of expected free energy (risk)
                %----------------------------------------------------------
                G(k) = G(k) + spm_kl_normal(m_q,c_q,m_p,c_p);

            end
        end

        % expected utility
        %------------------------------------------------------------------
        [~,j]  = max(F);
        W(:,4) = P{j}*W(:,4);

        % risk sensitive
        %------------------------------------------------------------------
        [~,j]  = min(G);
        W(:,5) = P{j}*W(:,5);


        % predictive posterior over RoR under this policy
        %------------------------------------------------------------------
        for t = 1:dT
            W_k    = W(:,5);
            ti     = i + t;
            Er(ti) = Ey(t,:)*W_k;
            Cr(ti) = W_k'*Cy{t}*W_k;
        end

    end

end

% legend labels
%--------------------------------------------------------------------------
spm_figure('GetWin','Portfolio management'); clf
try
    EFT = DEM.EFT;
catch
    for i = 1:m
        EFT{i} = sprintf('Asset %i',i);
    end
end

% Cumulative returns of policies
%--------------------------------------------------------------------------
xLim = [Ti - 1,Tf];
subplot(3,1,1), hold off
plot(exp(R)*100)
title('Cumulative returns of policies (%)','FontSize',14)
xlabel('time'), ylabel('%')
set(gca,'XLim',xLim)
legend(policies,'Location','northwest')

% Log returns of assets
%--------------------------------------------------------------------------
subplot(6,1,3), hold off
plot(L(:,1:m)*Ann)
title(sprintf('Rate of log return of assets (annualised)'),'FontSize',14)
ylabel('%'), set(gca,'XLim',xLim), legend(EFT)

% Allocation
%--------------------------------------------------------------------------
subplot(6,1,4), hold off
imagesc(Wpro)
title('Allocation: prospective','FontSize',14)
set(gca,'XLim',xLim)
set(gca,'YTickLabel', EFT)
set(gca,'YTick',1:m)

subplot(6,1,5), hold off
imagesc(Wret)
title('Allocation: retrospective','FontSize',14)
set(gca,'XLim',xLim)
set(gca,'YTickLabel', EFT)
set(gca,'YTick',1:m)

% Predicted and realised returns
%--------------------------------------------------------------------------
subplot(6,1,6), hold off
set(gca,'ColorOrderIndex',1)
spm_plot_ci(Er*Ann,Ann*Cr*Ann), hold on
set(gca,'ColorOrderIndex',1)
plot(D(:,5)*Ann,'.r')
title('Predicted and realised RoR','FontSize',14)
ylabel('%'), set(gca,'XLim',xLim)


% Annualised performance (D)
%==========================================================================
D     = exp(D(Ti:Tf,:)) - 1;
for i = 1:size(D,2)

    for t = 1:(size(D,1) - dT)
        ti    = t + (1:dT);

        % Annualised return
        %----------------------------------------------------------------------
        m(t)  = mean(D(ti,i))*Ann;

        % Volatility
        %----------------------------------------------------------------------
        s(t)  = std(D(ti,i))*Ann;

        % Sharpe ratio
        %----------------------------------------------------------------------
        r(t) = (m(t) - 2)/s(t);

    end

    % Annual return
    %----------------------------------------------------------------------
    tab{i,1} = mean(m);

    % Volatility
    %----------------------------------------------------------------------
    tab{i,2} = mean(s);

    % Sharpe ratio
    %----------------------------------------------------------------------
    tab{i,3} = mean(r);

    % Drawdown
    %----------------------------------------------------------------------
    tab{i,4} = min(D(:,i))*100;


end

VariableNames{1} = 'Annual RoR (%)';
VariableNames{2} = 'Volatility (%)';
VariableNames{3} = 'Sharpe ratio';
VariableNames{4} = 'Drawdown (%)';
RowNames = {'hold','max-ret','risk-ret','max-pro','risk-pro'};

Tab = cell2table(tab);
Tab.Properties.VariableNames = VariableNames;
Tab.Properties.RowNames      = RowNames

return

function [Ey,Cy,M] = spm_forecast_update(L,I,M,DEM,dT)
% FORMAT [Ey,Cy,M] = spm_forecast_update(L,I,M,DEM,dT)
% updates model parameters and generates forecast
% L   - rate of log return (per dt)
% U   - indicator variables (normalised)
% M   - flow model
% DEM - generative model
% dT  - forecast period
%__________________________________________________________________________

% get scale
%--------------------------------------------------------------------------
scale = diag(DEM.M(1).pE.scale);

% response and explanatory variables
%--------------------------------------------------------------------------
T    = size(L,1);               % current time
D    = dT;                      % past duration
t    = (T - D):T;               % past period
Y    = [I(t,:),L(t,:)];         % training data
X    = Y/scale;                 % scale

% model inversion with Variational Laplace
%==========================================================================
if false

    fig  = gcf;                     % get current figure
    M.X  = X;                       % legacy points in state-space
    B    = spm_ness_U(M);           % get state space and flow
    F    = gradient(M.X')';         % target flow

    % posterior over parameters
    %----------------------------------------------------------------------
    [Ep,Cp,Eh] = spm_nlsi_GN(M,B,F);

    % posterior precision of random fluctuations
    %----------------------------------------------------------------------
    W    = diag(exp(Eh));
    Ep.W = W;                       % posterior volatility
    M.W  = W;                       % update for flow model

    % Bayesian belief updating
    %----------------------------------------------------------------------
    M.pE = Ep;                      % and parameters of flow
    M.pC = Cp;                      % and covariance of parameters
    M.hE = Eh;                      % and log-precision

    % update initial expectations (parameters)
    %----------------------------------------------------------------------
    M.P  = Ep;

    % update generative model
    %======================================================================

    % model
    %----------------------------------------------------------------------
    DEM.M(1).pE = Ep;               % posterior esimates from VL
    DEM.M(1).pC = Cp;               % posterior esimates from VL
    DEM.M(1).V  = exp(16);
    DEM.M(1).W  = W;
    figure(fig)

end

% forecast
%==========================================================================
DEM.T      = T;
DEM.Y      = Y';
DEM.M(1).x = X(end,:)';
[Ey,Cy]    = spm_NESS_forecast(DEM,dT);

return

function [X] = spm_exp_conv(X,t)
% Epponential convolution
% FORMAT [X] = spm_exp_conv(X,t)
% X    - matrix
% sx   - kernel length
%__________________________________________________________________________
%
% spm_exp_conv is a onedimensional convolution of a matrix variable in
% working memory.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1999-2022 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
n  = size(X,1);
K  = exp(-(0:(6*t)).^2/(2*t^2));
C  = spm_convmtx(K(:),n);
C  = C(1:n,1:n);
C  = diag(1./sum(C,2))*C;
X  = C*X;






