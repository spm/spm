function Tab = spm_portfolio(L,I,DEM,Ti,dT)
% FORMAT spm_portfolio(L,DEM)
% simulate portfolio managment
%--------------------------------------------------------------------------
% L   - rate of log return per week
% I   - indicator variables
% DEM - generative model
% Ti  - time of initial investment
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
Tf    = size(L,1);              % time of final investment
m     = size(L,2);              % number of (aggregated) assets

% prior preferences
%--------------------------------------------------------------------------
try  RoR    = DEM.RoR;    catch, RoR    = 80;  end
try  Sharpe = DEM.Sharpe; catch, Sharpe = 1/2; end

m_p   = RoR/Ann;                % prior RoR per dt
s_p   = m_p/Sharpe;             % prior volitility (s.d.)
c_p   = s_p^2;                  % prior volitility (variance)

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
% W(1,:) = 0.1;  % 'VEA': 0.1
% W(2,:) = 0.05; % 'BIL': 0.05
% W(3,:) = 0.35; % 'AGG': 0.35
% W(4,:) = 0.05; % 'DBC': 0.05
% W(5,:) = 0.10; % 'VNQ': 0.10
% W(6,:) = 0.35; % 'SPY': 0.35

% simulate allocations and ensuing resturns
%==========================================================================
scale = diag(DEM.G.pE.scale);     % scaling of variables [I,L] -> X
TOL   = 128;                      % free energy threshold for rebalancing
R     = zeros(Tf,Np);             % cumulative RoR
S     = zeros(Tf,1);              % Surprise (Free energy)
D     = zeros(Tf,Np);             % RoR per dt
l     = (1:m) + size(I,2);        % indices of RoR
r     = -128:128;                 % domain of RoR (%)
for i = Ti:Tf

   
   % cumulative return on investment per dt
   %-----------------------------------------------------------------------
   D(i,:) = L(i,1:m)*W;
   R(i,:) = R(i - 1,:) + D(i,:);

   % Free energy
   %=======================================================================
   j   = (i - 2):i;               % immediate past
   X   = [I(j,:),L(j,:)]/scale;   % scaled states

   % place in generative model
   %-----------------------------------------------------------------------
   G   = DEM.G; G.Nmax = 1;       % generative model    
   G.x = X(1,:)';                 % initial state
   G.X = X;                       % legacy points in state-space

   % model inversion with Variational Laplace
   %-----------------------------------------------------------------------
   B   = spm_ness_U(G);           % basis functions
   f   = gradient(G.X')';         % flow

   % evaluate free energy
   %-----------------------------------------------------------------------
   [~,~,~,Fi] = spm_nlsi_GN(G,B,f);
   S(i)       = -Fi;


   % portfolio allocations
   %-----------------------------------------------------------------------
   W4(:,i)  = W(:,4);
   W5(:,i)  = W(:,5);

    % transactions
    %----------------------------------------------------------------------
    if ~rem(i,dT) || S(i) > TOL || i < (Ti - dT)

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
        [~,k]  = max(F);
        W(:,2) = P{k}*W(:,2);

        % risk sensitive
        %------------------------------------------------------------------
        [~,k]  = min(G);
        W(:,3) = P{k}*W(:,3);


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
            j           = 1:i;
            [Ez,Cz,DEM] = spm_forecast_update(L(j,:),I(j,:),DEM,dT + 1);
            
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

                % predictive posterior over outcomes
                %----------------------------------------------------------
                m_q  = Ey(t,:)*W_k;
                c_q  = W_k'*Cy{t}*W_k;

                % path integral of expected free energy (risk)
                %----------------------------------------------------------
                G(k) = G(k) + spm_kl_normal(m_q,c_q,m_p,c_p);

                % and drawdown contraint
                %----------------------------------------------------------
                f    = spm_Npdf(r,m_q*Ann,Ann*c_q*Ann);
                fc = 1;
                kl   = f*(log(f) - fc);
                G(k) = G(k) + kl;

            end
        end

        % expected utility
        %------------------------------------------------------------------
        [~,k]  = max(F);
        W(:,4) = P{k}*W(:,4);

        % risk sensitive
        %------------------------------------------------------------------
        [~,k]  = min(G);
        W(:,5) = P{k}*W(:,5);

        % predictive posterior over RoR under Risk-sensitive policy
        %------------------------------------------------------------------
        subplot(6,1,6), hold off

        % prior preferences
        %------------------------------------------------------------------
        f     = spm_Npdf(r,m_p*Ann,Ann*c_p*Ann);
        plot(r,f,'r'), hold on

        W_k   = W(:,5);
        for t = 1:dT

            % predictive posterior over outcomes
            %--------------------------------------------------------------
            m_q    = Ey(t,:)*W_k;
            c_q    = W_k'*Cy{t}*W_k;
            f      = spm_Npdf(r,m_q*Ann,Ann*c_q*Ann);
            plot(r,f,'b'), hold on

            ti     = i + t;
            Er(ti) = m_q;
            Cr(ti) = c_q;
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
title(sprintf('Rate of log return of assets (annualised)'),'FontSize',14)
ylabel('%'), set(gca,'XLim',xLim), legend(EFT)

% Allocation
%--------------------------------------------------------------------------
subplot(8,1,4), hold off
bar(W4(:,Ti:end)',1,'stacked','Edgecolor','none')
title('Expected Utility','FontSize',14)
legend(EFT)

subplot(8,1,5), hold off
bar(W5(:,Ti:end)',1,'stacked','Edgecolor','none')
title('Risk-sensitive','FontSize',14)
legend(EFT)

% Predicted and realised returns
%--------------------------------------------------------------------------
subplot(8,1,6), hold off
set(gca,'ColorOrderIndex',1)
spm_plot_ci(Er*Ann,Ann*Cr*Ann), hold on
set(gca,'ColorOrderIndex',1)
plot(D(:,5)*Ann,'.r')
title('Predicted and realised RoR','FontSize',14)
ylabel('%'), set(gca,'XLim',xLim)

% Free energy
%--------------------------------------------------------------------------
subplot(8,1,7), hold off
plot(S), hold on
plot([Ti,Tf],[1,1]*TOL,':r')
title('Uncertainty (Free Energy)','FontSize',14)
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

function [Ey,Cy,DEM] = spm_forecast_update(L,I,DEM,dT)
% FORMAT [Ey,Cy,DEM] = spm_forecast_update(L,I,DEM,dT)
% updates model parameters and generates forecast
% L   - rate of log return (per dt)
% U   - indicator variables (normalised)
% DEM - generative model
% dT  - forecast period
%__________________________________________________________________________

% get scale
%--------------------------------------------------------------------------
scale = diag(DEM.G.pE.scale);

% response and explanatory variables
%--------------------------------------------------------------------------
T    = size(L,1);                   % current time
D    = dT;                          % past duration
t    = (T - D):T;                   % past period
Y    = [I(t,:),L(t,:)];             % training data
X    = Y/scale;                     % scale

% model inversion with Variational Laplace
%==========================================================================
if true
    
    M    = DEM.G;                   % generative model of flow
    M.X  = X;                       % legacy points in state-space
    B    = spm_ness_U(M);           % get state space and flow
    F    = gradient(M.X')';         % target flow

    % posterior over parameters
    %----------------------------------------------------------------------
    [Ep,Cp] = spm_nlsi_GN(M,B,F);

    % Bayesian belief updating
    %----------------------------------------------------------------------
    DEM.G.pE = Ep;                  % parameters of flow
    DEM.G.pC = Cp;                  % and covariance of parameters

    % update model parameters in DEM
    %----------------------------------------------------------------------
    DEM.M(1).pE = Ep;               % posterior esimates from VL
    DEM.M(1).pC = Cp;               % posterior esimates from VL

end

% forecast
%==========================================================================
spm_figure('GetWin','Portfolio forecasts');
DEM.T      = T;
DEM.Y      = Y';
DEM.M(1).x = X(end,:)';
[Ey,Cy]    = spm_NESS_forecast(DEM,dT);

return
