function Tab = spm_portfolio(L,I,U,DEM,T,dT,Sr)
% FORMAT spm_portfolio(L,DEM)
% simulate portfolio managment
%--------------------------------------------------------------------------
% L   - rate of log return per week
% I   - indicator variables
% U   - exogenous variables
% DEM - generative model
% T   - time of inital investment
% dT  - time between rebalancing
% Sr  - Shaarpe ratio
%__________________________________________________________________________


% constraints on Bayesian belief updating
%--------------------------------------------------------------------------
DEM.G.pC = DEM.G.pC/16;

% default Sharpe ratio
%--------------------------------------------------------------------------
if nargin < 7, Sr = 3; end

% policies
%--------------------------------------------------------------------------
policies = {'Buy and hold', ...
    'Ex ante: expected RoR'...
    'Ex ante: risk sensitive'...
    'Ex post: Eexpected RoR'...
    'Ex post: risk sensitive'...
    };

% get sizes and times
%--------------------------------------------------------------------------
Ann   = 100*365/7;                % scaling for annualised % RoR
Ti    = T;                        % time of inital investment
Tf    = size(L,1);                % time of final investment
m     = size(L,2);                % number of (aggregated) assets
n     = size(I,2);                % number of indicator variables

% inital (flat) weights (i.e., ealth allocation)
%--------------------------------------------------------------------------
Np       = 5;                     % number of policies to simulate
W        = ones(6,Np,dT)/m;       % inital (flat) weights
% W(1,:,:) = 0.1;  % 'VEA': 0.1
% W(2,:,:) = 0.05; % 'BIL': 0.05
% W(3,:,:) = 0.35; % 'AGG': 0.35
% W(4,:,:) = 0.05; % 'DBC': 0.05
% W(5,:,:) = 0.10; % 'VNQ': 0.10
% W(6,:,:) = 0.35; % 'SPY': 0.35
W        = W(1:m,:,:);

% prior preferences
%--------------------------------------------------------------------------
Rr    = 256/Ann;                  % prior RoR percent
ni    = numel(Rr);                % prior preferences (Rate of return)
nj    = numel(Sr);                % prior preferences (Sharpe ratio)
m_p   = zeros(ni,nj);
c_p   = zeros(ni,nj);
for i = 1:ni
    for j = 1:nj
        m_p(i,j) = Rr(i);
        c_p(i,j) = (Rr(i)/Sr(j))^2;
    end
end


% Expected free energy: risk
%--------------------------------------------------------------------------
EFE = @(m_q,c_q,m_p,c_p) spm_kl_normal(m_q,c_q,m_p,c_p);

% Expected cost
%--------------------------------------------------------------------------
EC  = @(m_q) -m_q;

% prior preferences
%--------------------------------------------------------------------------
r   = -256:256;                   % domain of RoR (%) for plotting
% for i = 1:numel(m_p)
%     f = spm_Npdf(r,m_p(i)*Ann,Ann*c_p(i)*Ann);
%     plot(r,f,'b'), hold on
%     plot([0,0],    get(gca,'YLim'),'--r')
%     plot([-16,-16],get(gca,'YLim'),'--k')
% end

% graphics
%--------------------------------------------------------------------------
if ~isfield(DEM,'nograph')
    di  = datestr(DEM.G.date(Ti),'mmm-yy');
    df  = datestr(DEM.G.date(Tf),'mmm-yy');
    str = sprintf('Portfolio forecasts: %s to %s',di,df);
    spm_figure('GetWin',str); clf
end

% set up policies: moving funds from one assets to another nP times
%--------------------------------------------------------------------------
dP    = .5;                       % proportional transaction
cP    = .001;                     % proportional cost of transaction
nP    = 3;                        % number of transactions
p     = spm_combinations([m,m]);
np    = size(p,1);
Pk    = cell(1,np);
for k = 1:np
    i    = p(k,1);
    j    = p(k,2);
    Pij  = eye(m,m);

    % transfer (from asset j to i)
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
l     = (1:m) + n;                % indices of RoR
R     = zeros(Tf,Np,dT);          % cumulative RoR
S     = zeros(Tf,1);              % Surprise (Free energy)
D     = zeros(Tf,Np,dT);          % RoR per dt
E     = zeros(Tf,dT);             % expected free energy
for t = Ti:Tf

    % cumulative return on investment per dt
    %======================================================================
    fprintf('working on week %i\n',t - Ti + 1);
    DEM.T = t;
    for u = 1:dT
        D(t,:,u) = L(t,:)*W(:,:,u);
        R(t,:,u) = R(t - 1,:,u) + D(t,:,u);

        % plot outcomes
        %------------------------------------------------------------------
        if ~isfield(DEM,'nograph')
            subplot(3,2,6), hold on
            plot([1,1]*D(t,5,u)*Ann,get(gca,'YLim'),'g')
            plot([1,1]*R(t,5,u)*100,get(gca,'YLim'),'m')
            drawnow
        end

        % portfolio allocations
        %------------------------------------------------------------------
        W4(:,t,u)  = W(:,4,u);
        W5(:,t,u)  = W(:,5,u);

    end

    % Surpisal under this model
    %======================================================================
    % DEM  = spm_model_update(L,I,U,DEM,4);
    Ep   = DEM.M(1).pE;
    x    = DEM.M(1).gy([I(t,:),L(t,:)],U(t,:),Ep);
    S(t) = DEM.M(1).f(x,U(t,:),Ep,[],'S');

    % transaction stream
    %======================================================================
    w     = rem(t,dT) + 1;

    % if surpisal is large, rebalance all streams (or not)
    %----------------------------------------------------------------------
    if S(t) > 16
        %%% w = 1:dT;
    end

    for u = w

        % buy and hold
        %------------------------------------------------------------------
        W(:,1,u) = P{1}*W(:,1,u);

        % retrospective
        %==================================================================

        % exponential average: expected RoR
        %------------------------------------------------------------------
        % v   = DEM.M(1).E.v;
        v     = -1/4;
        a     = 1 - exp(v);
        EL    = a*L(t,:);
        for i = 1:32
            EL = EL + a*(1 - a)^i*L(t - i,:);    % expected RoR
        end
        CL    = cov(L((t - 32):t,:));            % covariance of RoR
        F     = zeros(1,np);
        G     = zeros(ni,nj,np);
        for k = 1:np

            % expected utility
            %--------------------------------------------------------------
            W_k  = P{k}*W(:,2,u);
            m_q  = EL*W_k;

            % path integral of expected free energy (risk)
            %--------------------------------------------------------------
            F(k) = F(k) - m_q;

            % risk sensitive
            %--------------------------------------------------------------
            W_k  = P{k}*W(:,3,u);
            m_q  = EL*W_k;
            c_q  = W_k'*CL*W_k;

            % for every prior preference
            %--------------------------------------------------------------
            for i = 1:ni
                for j = 1:nj
                    kl       = EFE(m_q,c_q,m_p(i,j),c_p(i,j));
                    G(i,j,k) = G(i,j,k) + kl;
                end
            end

        end

        % expected utility
        %------------------------------------------------------------------
        [~,k]    = min(F);
        W(:,2,u) = P{k}*W(:,2,u);

        % risk sensitive
        %------------------------------------------------------------------
        [~,j]    = min(G,[],'all');
        [~,~,k]  = ind2sub([ni,nj,np],j);
        W(:,3,u) = P{k}*W(:,3,u);


        % prospective
        %==================================================================
        if false

            % upper bound on performance, if the future were known
            %--------------------------------------------------------------
            nT    = dT;
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
            [Ez,Cz,DEM] = spm_forecast_update(L(s,:),I(s,:),U,DEM,dT);

            nT    = numel(Cz) - 1;
            for s = 1:nT
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
            W_k   = P{k}*W(:,4,u);
            for s = 1:nT

                % predictive posterior over outcomes
                %----------------------------------------------------------
                m_q  = Ey(s,:)*W_k;

                % path integral of expected free energy (risk)
                %----------------------------------------------------------
                F(k) = F(k) - m_q;

            end

            % risk sensitive
            %--------------------------------------------------------------
            W_k   = P{k}*W(:,5,u);
            for s = 1:nT

                % predictive posterior over outcomes
                %----------------------------------------------------------
                m_q  = Ey(s,:)*W_k;
                c_q  = W_k'*Cy{s}*W_k;

                % path integral of expected free energy (risk)
                %----------------------------------------------------------
                for i = 1:ni
                    for j = 1:nj
                        kl       = EFE(m_q,c_q,m_p(i,j),c_p(i,j));
                        G(i,j,k) = G(i,j,k) + kl;
                    end
                end
            end

        end

        % expected utility
        %------------------------------------------------------------------
        [~,k]    = min(F);
        W(:,4,u) = P{k}*W(:,4,u);

        % risk sensitive: find best policy (k) under preferences (j)
        %------------------------------------------------------------------
        [e,j]    = min(G,[],'all');
        [i,j,k]  = ind2sub([ni,nj,np],j);
        K        = P{k};

        % or Bayesian model averaging
        %------------------------------------------------------------------
        % pk    = -G*16;
        % pk(:) = spm_softmax(pk(:));
        % pk    = sum(pk,[1,2]);
        % K     = 0;
        % for k = 1:np
        %     K = K + P{k}*pk(k);
        % end
        %------------------------------------------------------------------
        W(:,5,u) = K*W(:,5,u);

        % predictive posterior over RoR under Risk-sensitive policy
        %==================================================================
        subplot(3,2,6), hold off

        % prior preferences (j)
        %------------------------------------------------------------------
        if ~isfield(DEM,'nograph')
            f = spm_Npdf(r,m_p(i,j)*Ann,Ann*c_p(i,j)*Ann);
            plot(r,f,'r'), hold on
            plot([0,0],get(gca,'YLim'),'--r')
            axis square
        end

        W_k   = W(:,5,u);
        for s = 1:nT

            % predictive posterior over outcomes
            %--------------------------------------------------------------
            m_q  = Ey(s,:)*W_k;
            c_q  = W_k'*Cy{s}*W_k;
            f    = spm_Npdf(r,m_q*Ann,Ann*c_q*Ann);

            if ~isfield(DEM,'nograph')
                if s > 1
                    plot(r,f,'c'), hold on
                else
                    plot(r,f,'b'), hold on
                end
                title('Predictive and preferred densities')
                xlabel('Annualised RoR (%)')
            end

            % predictive posterior over RoR
            %--------------------------------------------------------------
            Er(t + s,u) = m_q;
            Cr(t + s,u) = c_q;
        end

        % expected free energy
        %------------------------------------------------------------------
        E(t,u) = e;
    end

end

% average over phase (portfolio streams)
%--------------------------------------------------------------------------
R   = mean(R,3);
D   = mean(D,3);
W4  = mean(W4,3);
W5  = mean(W5,3);
Er  = mean(Er,2);
Cr  = mean(Cr,2);

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

% Predicted and realised returns
%--------------------------------------------------------------------------
subplot(8,1,6), hold off
set(gca,'ColorOrderIndex',1)
spm_plot_ci(Er'*Ann,Ann*Cr'*Ann), hold on
set(gca,'ColorOrderIndex',1)
plot(D(:,5)*Ann,'.b','MarkerSize',8), hold on
plot(xLim,[0,0],':b')
title('Predicted and realised RoR','FontSize',12)
ylabel('%'), set(gca,'XLim',xLim)

% Variational free energy
%--------------------------------------------------------------------------
subplot(8,1,7), hold off
plot(S),  hold on
plot([Ti,Tf],[1,1]*5,':r')
title('Uncertainty: Surprisal','FontSize',12)
ylabel('nats'), set(gca,'XLim',xLim)

% Expected free energy
%--------------------------------------------------------------------------
subplot(8,1,8), hold off
plot(E,'.','MarkerSize',16), hold on
plot([Ti,Tf],[1,1]*3,':r')
title('Confidence: Expected Free Energy','FontSize',12)
ylabel('nats'), set(gca,'XLim',xLim)


% Annualised performance (D)
%==========================================================================
D     = exp(D(Ti:Tf,:)) - 1;
for i = 1:size(D,2)

    % monthly annualised return
    %----------------------------------------------------------------------
    for t = 1:(size(D,1) - 4)
        ti    = t + (1:4);
        m(t)  = mean(D(ti,i))*Ann;
    end

    % Annual return
    %----------------------------------------------------------------------
    tab(i,1) = mean(D(:,i))*Ann;

    % Volatility
    %----------------------------------------------------------------------
    tab(i,2) = std(m);

    % Sharpe ratio
    %----------------------------------------------------------------------
    tab(i,3) = (mean(m) - 2)/std(m);

    % Drawdown
    %----------------------------------------------------------------------
    tab(i,4) = min(D(:,i))*100;

end

VariableNames{1} = 'Annual RoR (%)';
VariableNames{2} = 'Volatility (%)';
VariableNames{3} = 'Sharpe ratio';
VariableNames{4} = 'Drawdown (%)';
RowNames = {'hold','ex post EV','ex post KL','ex ante EV','ex ante KL'};

Tab = array2table(tab);
Tab.Properties.VariableNames = VariableNames;
Tab.Properties.RowNames      = RowNames;

return


function [Ey,Cy,DEM] = spm_forecast_update(L,I,U,DEM,dT)
% FORMAT [Ey,Cy,DEM] = spm_forecast_update(L,I,U,DEM,dT)
% generates forecast
% L   - rate of log return (per dt)
% I   - indicator variables
% U   - exogenous variables
% DEM - generative model
% dT  - forecast period
%__________________________________________________________________________

% get parameters
%--------------------------------------------------------------------------
pE   = DEM.G.pE;
gy   = DEM.M(1).gy;

% response and explanatory variables
%--------------------------------------------------------------------------
T    = DEM.T;                          % current time
x    = gy([I(T,:),L(T,:)],U(T,:),pE);  % recent states
u    = U(T,:);

% forecast
%==========================================================================
t          = T:(T + dT);
DEM.U      = U(t,:)';                  % future causes
DEM.Y      = [I(1:T,:),L(1:T,:)]';     % past data
DEM.M(1).x = x(:);                     % intial state
DEM.M(2).v = u(:);                     % intial cause
[Ey,Cy]    = spm_NESS_forecast(DEM);

return

function [DEM] = spm_model_update(L,I,U,DEM,dT)
% FORMAT [DEM] = spm_model_update(L,I,U,DEM,dT)
% updates model parameters
% L   - rate of log return (per dt)
% I   - indicator variables
% U   - exogenous variables
% DEM - generative model
% dT  - forecast period
%__________________________________________________________________________

% get parameters and observer function
%--------------------------------------------------------------------------
pE   = DEM.G.pE;
gy   = DEM.M(1).gy;

% response and explanatory variables
%--------------------------------------------------------------------------
T    = DEM.T;                       % current time
t    = (T - dT):T;                  % recent past
Y    = [I(t,:),L(t,:)];             % recent data
X    = gy(Y,U(t,:),pE);             % recent states

% model inversion with Variational Laplace
%==========================================================================
M    = DEM.G;                       % generative model of flow
M.X  = X;                           % legacy points in state-space
f    = gradient(X')';               % target flow
B    = spm_ness_U(M);               % get basis functions

% posterior over parameters
%--------------------------------------------------------------------------
M.nograph = 1;
M.Nmax    = 4;
[Ep,Cp]   = spm_nlsi_GN(M,B,f);

% Bayesian belief updating
%--------------------------------------------------------------------------
M.pE = Ep;                          % parameters of flow
M.pC = Cp;                          % and covariance of parameters

% update model parameters in DEM
%--------------------------------------------------------------------------
DEM.M(1).pE = Ep;                   % posterior esimates from VL
DEM.M(1).pC = Cp;                   % posterior esimates from VL
DEM.G       = M;                    % generative model of flow

return

% Alternative update based upon generalised filtering
%==========================================================================

% response and explanatory variables
%--------------------------------------------------------------------------
T    = DEM.T;                       % current time
t    = (T - dT):T;                  % recent past
Y    = [I(t,:),L(t,:)];             % recent data
u    = [U(t,:)];                    % exogenous inputs

DEM.Y = Y';
DEM.U = u';
DEM   = rmfield(DEM,'X');


% get scale
%--------------------------------------------------------------------------
pC          = spm_zeros(DEM.M(1).pE);
pC.trend(:) = ones(size(pC.trend))/32;
pC.W        = eye(size(pC.W));
DEM.M(1).pC = pC;

DEM         = spm_LAP(DEM);

DEM.M(1).pE = DEM.M(1).Ep;
DEM.G.pE    = DEM.M(1).Ep;

return

% Alternative update based updating scale and trend parameters
%==========================================================================

% response and explanatory variables
%--------------------------------------------------------------------------
T    = DEM.T;                       % current time
t    = (T - dT):T;                  % recent past
Y    = [I(t,:),L(t,:)];             % recent data
u    = [U(t,:)];                    % exogenous inputs

% get scale
%--------------------------------------------------------------------------
P       = DEM.M(1).pE;
P.trend = u\Y;
P.scale = sqrtm(cov(Y - u*P.trend));

DEM.M(1).pE = P;
DEM.G.pE    = P;

return