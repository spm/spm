function dSYdP = spm_COVID_R_cii(DCM,U,name)
% Graphics for coronavirus simulations - with confidence intervals
% FORMAT dSYdP = spm_COVID_R_ci(DCM,U)
% DCM.Ep - posterior expectations
% DCM.Cp - posterior covariances
% DCM.Y  - empirical data
% DCM.M  - model
%
% U      - output to evaluate [default: 1]
%
% dSYdP  - first-order sensitivity (with respect to outcome U)
%
% This routine evaluates a trajectory of outcome variables from a COVID
% model and plots the expected trajectory and accompanying Bayesian
% credible intervals (of 90%). If empirical data are supplied, these will
% be overlaid on the confidence intervals. By default, 365 days are
% evaluated. In addition, posterior and prior expectations are provided in
% a panel. this confidence interval potting routine handles multiple region
% models and returns both a sensitivity analysis and posterior predictive
% density over specified outcomes (in U).
%
% Although the covid model is non-linear in the parameters, one can use a
% first-order Taylor expansion to evaluate the confidence intervals in
% terms of how the outcomes change with parameters. This, in combination
% with the well-known overconfidence of variational inference, usually
% requires a slight inflation of uncertainty. Here, the posterior
% covariance is multiplied by a factor of four.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% setup and default: number of outcomes to evaluate
%--------------------------------------------------------------------------
if nargin < 1, U = 1;     end
if nargin < 2, name = ''; end

Ep  = DCM.Ep;                                      % posterior expectations
Cp  = DCM.Cp;                                      % posterior covariances
Z   = DCM.Y;                                       % empirical data
M   = DCM.M;                                       % model


% priors and names for plotting
%--------------------------------------------------------------------------
[pE,pC,out] = spm_COVID_priors;
[pE,pC,str] = spm_COVID_priors_R(DCM.M.data);

% times (and compensate for overconfidence)
%--------------------------------------------------------------------------
M.T = 365;
t   = (1:M.T)/7;

% evaluate confidence intervals (using a Taylor expansion)
%==========================================================================

% dimension reduction of parameter space
%--------------------------------------------------------------------------
V        = spm_svd(Cp,0);
[dYdV,Y] = spm_diff(@(P,M,U)spm_COVID_US(P,M,U),Ep,M,U,1,{V});

% conditional covariances for each country
%--------------------------------------------------------------------------
Nk    = size(Y,3);
C     = cell(Nk,1);
s     = zeros(Nk,1);
SS    = 0;
dSSdP = 0;
for k = 1:Nk
    
    % trajectories of death rates
    %----------------------------------------------------------------------
    try
        for j = 1:numel(dYdV)
            D{j} = dYdV{j}(:,1,k);
        end
        dYdP  = spm_cat(D)*V';
    catch
        dYdP  = dYdV*V';
    end
    C{k}  = diag(16*Y(:,1,k));
    
    % cumulative death rates per region
    %----------------------------------------------------------------------
    S     = cumsum(Y(:,1,k));
    
    % cumulative death rates over regions
    %----------------------------------------------------------------------
    SS    = SS + S;
    dSSdP = dSSdP + dYdP;
    s(k)  = S(end);
    
end
CSS  = 16*SS;
cs   = 16*s;

% Sensitivity analysis: which factors determine cumulative deaths?
%==========================================================================
spm_figure('GetWin',['Sensitivity' name]); clf
%--------------------------------------------------------------------------
% (sensitivity analysis) This figure reports the effect of changing each
% parameter on cumulative deaths over an 18-month period. The upper panel
% shows the rate of increase (or decrease) in cumulative deaths per unit
% change in the (log) parameters. These sensitivity metrics are based upon
% a first order Taylor expansion about the maximum a posteriori values
% shown in the lower panel. The blue bars correspond to the most likely
% estimate and the pink bars report the 90% credible intervals.
%--------------------------------------------------------------------------
ii    = ismember([str.names],'erc');
ep    = spm_vec(Ep);
dSYdP = sum(dSSdP(:,ii));
subplot(3,1,1)
bar(dSYdP)
xlabel('connectivity parameter')
ylabel('sensitivity')
title('First-order sensitivity','FontSize',16), box off

subplot(3,1,2), hold off
spm_plot_ci(ep(ii),Cp(ii,ii),[],[],'exp')
xlabel('connectivity parameter')
ylabel('probability of leaving a State')
title('MAP connectivity estimates','FontSize',16), box off


% projections and confidence intervals
%==========================================================================
spm_figure('GetWin',['Projections' name]); clf
%--------------------------------------------------------------------------
% (projected mortality rates) This figure reports the predicted trajectory
% of death rates for each state. The upper panel shows the predicted time
% course of death rates per day. The expected rates are shown as blue
% lines, while the shaded blue areas correspond to 90% Bayesian credible
% intervals. The dots report empirical data observed to date. The lower
% left panel shows the cumulative deaths based upon these predictions, in
% terms of the total expectation (blue bars) and accompanying 90% credible
% intervals (pink bars). The lower right panel shows the cumulative deaths
% over all regions, in terms of the posterior expectation (blue line) and
% confidence intervals (shaded area). The (black dots) correspond to
% empirical data. Note that variational procedures are notoriously
% overconfident. We have compensated for this here by increasing the
% confidence intervals by a factor of two; however, one should view these
% confidence intervals as overconfident in and of themselves.
%--------------------------------------------------------------------------
subplot(3,2,1), hold off
for k = 1:Nk
    subplot(3,2,1)
    [y,i] = max(Y(:,U,k));
    spm_plot_ci(Y(:,U,k)',C{k},t), hold on    
    text(t(i),y,str.regions{k},'FontSize',8,'Color','b')
    try, plot(t(1:numel(Z(:,U,k))),Z(:,U,k),'.k'), end
end
xlabel('time (weeks)'),ylabel('number of deaths per day')
title(out.outcome(U),'FontSize',16)
box off, spm_axis tight
set(gca,'YLim',[0 max(max(Y(:,U,:)))]);

subplot(3,2,2), hold off
for k = 1:Nk
    subplot(3,2,2)
    [y,i] = max(Y(:,U,k));
    spm_plot_ci(Y(:,U,k)',C{k},t), hold on
    try, plot(t(1:numel(Z(:,U,k))),Z(:,U,k),'.k'), end
end
xlabel('time (weeks)'),ylabel('number of deaths per day')
title(out.outcome(U),'FontSize',16)
box off, spm_axis tight
set(gca,'YLim',[0 4*min(max(Y(:,U,:)))]);

% cumulative deaths
%--------------------------------------------------------------------------
subplot(3,2,3), hold off
spm_plot_ci(s,cs)
ylabel('number of deaths')
set(gca,'XTick',1:numel(str.regions),...
        'Xticklabel',str.regions,...
        'XTickLabelRotation',90,'FontSize',10)
title('Cumulative deaths','FontSize',16)
axis square, box off, spm_axis tight

N     = [DCM.M.data.pop]'/100;
subplot(3,2,4), hold off
spm_plot_ci(s./N,cs./(N.^2)), hold on
plot([1 numel(s)],[1,1]/10,'-.')
ylabel('cumulative deaths per capita (%)')
set(gca,'XTick',1:numel(str.regions),...
        'Xticklabel',str.regions,...
        'XTickLabelRotation',90,'FontSize',10)
title('Mortality rate (%)','FontSize',16)
axis square, box off, spm_axis tight

subplot(3,1,3), hold off
spm_plot_ci(SS',CSS,t), hold on
try, plot(t(1:size(Z,1)),cumsum(sum(Z(:,U,:),3)),'.k'), end
xlabel('time (weeks)'),ylabel('number of cases')
title('Cumulative deaths','FontSize',16)
axis square, box off, spm_axis tight


return








