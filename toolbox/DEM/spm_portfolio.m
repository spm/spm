function Tab = spm_portfolio(L,I,U,DEM,dT)
% FORMAT spm_portfolio(L,DEM)
% simulate portfolio managment
%--------------------------------------------------------------------------
% L   - rate of log return per week
% I   - indicator variables
% U   - exogenous variables
% DEM - generative model
% dT  - (maxiumum) time between rebalancing
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
dt    = DEM.G.dt;               % period (days)
Ann   = 100*365/dt;             % scaling for annualised % RoR
Ti    = DEM.T;                  % time of inital investment
Tf    = size(L,1);              % time of final investment
m     = size(L,2);              % number of (aggregated) assets
DEM.t = Ti;

% prior preferences
%--------------------------------------------------------------------------
RoR   = [8, 16, 32, 64, 128, 256];    % prior RoR percent
Rat   = [2, 3, 4];                    % Sharpe ratio
ni    = numel(RoR);             % number of prior preferences
nj    = numel(Rat);             % number of prior preferences
m_p   = zeros(ni,nj);
c_p   = zeros(ni,nj);
p_p   = zeros(ni,nj);
for i = 1:ni
    for j = 1:nj
        m_p(i,j) = RoR(i)/Ann;
        c_p(i,j) = ((RoR(i) + 16)/Ann/Rat(j))^2;
        p_p(i,j) = -log(RoR(i)) - log(Rat(j));

    end
end
p_p   = p_p - min(p_p(:));
p_p   = p_p*0; %%%



% set up policies: moving funds from one assets to another nP times
%--------------------------------------------------------------------------
dP    = .5;                     % proportional transaction
cP    = .001;                   % proportional cost of transaction
nP    = 3;                      % number of transactions
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

% inital (flat) weights (i.e., ealth allocation)
%--------------------------------------------------------------------------
Np     = 5;                       % number of policies to simulate
W      = ones(m,Np)/m;            % inital (flat) weights
% or
W(1,:) = 0.1;  % 'VEA': 0.1
W(2,:) = 0.05; % 'BIL': 0.05
W(3,:) = 0.35; % 'AGG': 0.35
W(4,:) = 0.05; % 'DBC': 0.05
W(5,:) = 0.10; % 'VNQ': 0.10
W(6,:) = 0.35; % 'SPY': 0.35


% amplitude of random fluctuations
%--------------------------------------------------------------------------
l     = (1:m) + size(I,2);        % indices of RoR
scale = diag(DEM.G.pE.scale);     % scaling of variables [I,L] -> X

% simulate allocations and ensuing resturns
%==========================================================================
TOL   = 64;                       % free energy threshold for rebalancing
r     = -128:256;                 % domain of RoR (%) for plotting
R     = zeros(Tf,Np);             % cumulative RoR
S     = zeros(Tf,1);              % Surprise (Free energy)
D     = zeros(Tf,Np);             % RoR per dt
E     = [];                       % expected free energy
T     = [];                       % times of policy evaluation
for t = Ti:Tf

    % cumulative return on investment per dt
    %======================================================================
    D(t,:) = L(t,:)*W;
    R(t,:) = R(t - 1,:) + D(t,:);

    % plot outcomes
    %----------------------------------------------------------------------
    subplot(6,2,12), hold on
    plot([1,1]*D(t,5)*Ann,get(gca,'YLim'),'g')
    plot([1,1]*R(t,5)*100,get(gca,'YLim'),'m')


    % Free energy
    %======================================================================
    j   = (t - 2):t;               % immediate past
    X   = [I(j,:),L(j,:)]/scale;   % scaled states

    % place in generative model
    %----------------------------------------------------------------------
    G   = DEM.G; G.Nmax = 1;       % generative model
    G.x = X(1,:)';                 % initial state
    G.X = X;                       % legacy points in state-space

    % model inversion with Variational Laplace
    %----------------------------------------------------------------------
    B   = spm_ness_U(G);           % basis functions
    B.u = U(j,:);                  % exogenous input
    f   = gradient(G.X')';         % flow

    % evaluate free energy
    %----------------------------------------------------------------------
    [~,~,~,Fi] = spm_nlsi_GN(G,B,f);
    S(t)       = -Fi;

    % portfolio allocations
    %----------------------------------------------------------------------
    W4(:,t)  = W(:,4);
    W5(:,t)  = W(:,5);

    % transactions
    %======================================================================
    if ~rem(t,dT) || S(t) - S(t - 1) > TOL

        % buy and hold
        %------------------------------------------------------------------
        W(:,1)  = P{1}*W(:,1);

        % retrospective
        %==================================================================
        EL    = L(t,:);                     % expected RoR
        CL    = cov(L((t - 16):t,:));       % sample covriance of RoR
        F     = zeros(1,np);
        G     = zeros(ni,nj,np);
        for k = 1:np

            % expected utility
            %--------------------------------------------------------------
            W_k  = P{k}*W(:,2);
            m_q  = EL*W_k;

            % path integral of expected free energy (risk)
            %--------------------------------------------------------------
            kl   = m_q;
            F(k) = F(k) - kl;

            % risk sensitive
            %--------------------------------------------------------------
            W_k  = P{k}*W(:,3);
            m_q  = EL*W_k;
            c_q  = W_k'*CL*W_k;

            % for every prior preference
            %--------------------------------------------------------------
            for i = 1:ni
                for j = 1:nj
                    kl       = spm_kl_normal(m_q,c_q,m_p(i,j),c_p(i,j));
                    G(i,j,k) = G(i,j,k) + kl + p_p(i,j);
                end
            end

        end

        % expected utility
        %------------------------------------------------------------------
        [~,k]  = min(F);
        W(:,2) = P{k}*W(:,2);

        % risk sensitive
        %------------------------------------------------------------------
        [~,j]   = min(G,[],'all');
        [~,~,k] = ind2sub([ni,nj,np],j);
        W(:,3)  = P{k}*W(:,3);


        % prospective
        %==================================================================
        if false

            % upper bound on performance, if the future were known
            %--------------------------------------------------------------
            for s = 1:dT
                try
                    Ey(s,:) = L(s + t,1:m);
                    Cy{s}   = cov(L((1:dT) + s,:));
                catch
                    Ey(s,:) = L(end,1:m);
                    Cy{s}   = cov(L(end - (1:dT) + 1,:));
                end
            end

        else

            % posterior predictive density over RoR
            %--------------------------------------------------------------
            s           = 1:t;
            [Ez,Cz,DEM] = spm_forecast_update(L(s,:),I(s,:),U(s,:),DEM,dT + 1);

            for s = 1:dT
                Ey(s,:) = Ez(l,s + 1)';
                Cy{s}   = Cz{s + 1}(l,l);
            end

        end


        % for each policy
        %------------------------------------------------------------------
        F     = zeros(1,np);
        G     = zeros(ni,nj,np);
        for k = 1:np

            % expected utility
            %--------------------------------------------------------------
            W_k   = P{k}*W(:,4);
            for s = 1:dT

                % predictive posterior over outcomes
                %----------------------------------------------------------
                m_q  = Ey(s,:)*W_k;

                % path integral of expected free energy (risk)
                %----------------------------------------------------------
                kl   = m_q;
                F(k) = F(k) - kl;

            end

            % risk sensitive
            %--------------------------------------------------------------
            W_k   = P{k}*W(:,5);
            for s = 1:dT

                % predictive posterior over outcomes
                %----------------------------------------------------------
                m_q  = Ey(s,:)*W_k;
                c_q  = W_k'*Cy{s}*W_k;

                % path integral of expected free energy (risk)
                %----------------------------------------------------------
                for i = 1:ni
                    for j = 1:nj
                        kl       = spm_kl_normal(m_q,c_q,m_p(i,j),c_p(i,j));
                        G(i,j,k) = G(i,j,k) + kl + p_p(i,j);
                    end
                end
            end

        end

        % expected utility
        %------------------------------------------------------------------
        [~,k]  = min(F);
        W(:,4) = P{k}*W(:,4);

        % risk sensitive: find best policy (k) under preferences (j)
        %------------------------------------------------------------------
        [e,j]   = min(G,[],'all');
        [i,j,k] = ind2sub([ni,nj,np],j);
        W(:,5)  = P{k}*W(:,5);

        % predictive posterior over RoR under Risk-sensitive policy
        %==================================================================
        subplot(6,2,12), hold off

        % prior preferences (j)
        %------------------------------------------------------------------
        f     = spm_Npdf(r,m_p(i,j)*Ann,Ann*c_p(i,j)*Ann);
        plot(r,f,'r'), hold on
        plot([0,0],get(gca,'YLim'),'--r')

        W_k   = W(:,5);
        for s = 1:dT

            % predictive posterior over outcomes
            %--------------------------------------------------------------
            m_q  = Ey(s,:)*W_k;
            c_q  = W_k'*Cy{s}*W_k;
            f    = spm_Npdf(r,m_q*Ann,Ann*c_q*Ann);
            if s > 1
                plot(r,f,'b:'), hold on
            else
                plot(r,f,'b'), hold on
            end
            title('Predictive and preferred densities')
            xlabel('Annualised RoR (%)')

            Er(s + t) = m_q;
            Cr(s + t) = c_q;
        end

        subplot(6,4,21)
        imagesc(-p_p)
        xlabel('Sharpe'), ylabel('RoR'), title('Prior preferences')

        subplot(6,4,22)
        imagesc(-sum(G,3))
        xlabel('Sharpe'), ylabel('RoR'), title('Expected free energy')

        % expected free energy
        %------------------------------------------------------------------
        E(end + 1) = e;
        T(end + 1) = t;

    end

end

% legend labels
%--------------------------------------------------------------------------
di  = datestr(DEM.G.date(Ti),'mmm-yy');
df  = datestr(DEM.G.date(Tf),'mmm-yy');
str = sprintf('Portfolio management: %s to %s',di,df);
spm_figure('GetWin',str); clf
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
subplot(4,1,1), hold off
plot(exp(R)*100)
title('Cumulative returns of policies (%)','FontSize',14)
xlabel('time'), ylabel('%')
set(gca,'XLim',xLim)
legend(policies,'Location','northwest')

% Log returns of assets
%--------------------------------------------------------------------------
subplot(8,1,3), hold off
plot(L(:,1:m)*Ann)
title(sprintf('Rate of log return of assets (annualised)'),'FontSize',12)
ylabel('%'), set(gca,'XLim',xLim), legend(EFT)

% Allocation
%--------------------------------------------------------------------------
subplot(8,1,4), hold off
bar(W4(:,Ti:end)',1,'stacked','Edgecolor','none')
title('Expected Utility','FontSize',12)
legend(EFT)

subplot(8,1,5), hold off
bar(W5(:,Ti:end)',1,'stacked','Edgecolor','none')
title('Risk-sensitive','FontSize',12)
legend(EFT)

% Predicted and realised returns
%--------------------------------------------------------------------------
subplot(8,1,6), hold off
set(gca,'ColorOrderIndex',1)
spm_plot_ci(Er*Ann,Ann*Cr*Ann), hold on
set(gca,'ColorOrderIndex',1)
plot(D(:,5)*Ann,'.r')
title('Predicted and realised RoR','FontSize',12)
ylabel('%'), set(gca,'XLim',xLim)

% Variational free energy
%--------------------------------------------------------------------------
subplot(8,1,7), hold off
plot(S), hold on
plot([Ti,Tf],[1,1]*TOL,':r')
title('Uncertainty: Variational Free Energy','FontSize',12)
ylabel('nats'), set(gca,'XLim',xLim)

% Expected free energy
%--------------------------------------------------------------------------
subplot(8,1,8), hold off
plot(T,E,'o'), hold on
plot([Ti,Tf],[1,1]*3,':r')
title('Confidence: Expected Free Energy','FontSize',12)
ylabel('nats'), set(gca,'XLim',xLim)


% Annualised performance (D)
%==========================================================================
D     = exp(D(Ti:Tf,:)) - 1;
dt    = 4;
for i = 1:size(D,2)

    % monthly annualised return
    %----------------------------------------------------------------------
    for t = 1:(size(D,1) - dt)
        ti    = t + (1:dt);
        m(t)  = mean(D(ti,i))*Ann;
    end

    % Annual return
    %----------------------------------------------------------------------
    tab(i,1) = mean(m);

    % Volatility
    %----------------------------------------------------------------------
    tab(i,2) = std(m);

    % Sharpe ratio
    %----------------------------------------------------------------------
    tab(i,3) = (mean(m) - 1)/std(m);

    % Drawdown
    %----------------------------------------------------------------------
    tab(i,4) = min(D(:,i))*100;

end

VariableNames{1} = 'Annual RoR (%)';
VariableNames{2} = 'Volatility (%)';
VariableNames{3} = 'Sharpe ratio';
VariableNames{4} = 'Drawdown (%)';
RowNames = {'hold','max-ret','risk-ret','max-pro','risk-pro'};

Tab = array2table(tab);
Tab.Properties.VariableNames = VariableNames;
Tab.Properties.RowNames      = RowNames;

return

function [Ey,Cy,DEM] = spm_forecast_update(L,I,U,DEM,dT)
% FORMAT [Ey,Cy,DEM] = spm_forecast_update(L,I,U,DEM,dT)
% updates model parameters and generates forecast
% L   - rate of log return (per dt)
% I   - indicator variables
% U   - exogenous variables
% DEM - generative model
% dT  - forecast period
%__________________________________________________________________________

% get scale
%--------------------------------------------------------------------------
scale = diag(DEM.G.pE.scale);

% response and explanatory variables
%--------------------------------------------------------------------------
T    = size(L,1);                   % current time
t    = (T - dT):T;                  % recent past
Y    = [I(t,:),L(t,:)];             % recent data
X    = Y/scale;                     % recent states

% model inversion with Variational Laplace
%==========================================================================
if true

    M    = DEM.G;                   % generative model of flow
    M.X  = X;                       % legacy points in state-space
    B    = spm_ness_U(M);           % get state space and flow
    B.u  = U(t,:);
    F    = gradient(M.X')';         % target flow

    % posterior over parameters
    %----------------------------------------------------------------------
    M.nograph = 1;
    M.Nmax    = 32;
    [Ep,Cp]   = spm_nlsi_GN(M,B,F);

    % Bayesian belief updating
    %----------------------------------------------------------------------
    M.pE = Ep;                      % parameters of flow
    M.pC = Cp;                      % and covariance of parameters

    % update model parameters in DEM
    %----------------------------------------------------------------------
    DEM.M(1).pE = Ep;               % posterior esimates from VL
    DEM.M(1).pC = Cp;               % posterior esimates from VL
    DEM.G       = M;                % generative model of flow

end

% forecast
%==========================================================================
df  = datestr(DEM.G.date(DEM.t),'mmm-yy');
str = sprintf('Portfolio forecasts: from %s',df);
spm_figure('GetWin',str);

DEM.Y      = Y';                    % update recent data
DEM.T      = T;                     % update current time
DEM.M(1).x = X(end,:)';             % intial state
DEM.M(2).v = U(end,:)';             % intial cause

[Ey,Cy]    = spm_NESS_forecast(DEM,dT);

return
