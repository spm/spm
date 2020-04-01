function spm_COVID_ci(Ep,Cp,Z)
% Graphics for coronavirus simulations - with confidence intervals
% FORMAT spm_COVID_plot(Y,X,Z)
% Ep     - posterior expectations
% Cp     - posterior covariances
% Z      - optional empirical data
%
% This routine evaluates a trajectory of outcome variables from a COVID
% model and plots the expected trajectory and accompanying Bayesian
% credible intervals (of 90%). If empirical data are supplied, these will
% be overlaid on the confidence intervals. By default, hundred and 28 days
% are evaluated. In addition, posterior and prior expectations are provided
% in a panel.
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
% $Id: spm_COVID_ci.m 7810 2020-04-01 13:58:56Z spm $

% priors and names for plotting
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% compensate for (variational) overconfidence
%--------------------------------------------------------------------------
Cp      = Cp*4;

% evaluate confidence intervals (using a Taylor expansion)
%==========================================================================

% changes in outcomes with respect to parameters
%--------------------------------------------------------------------------
M.T      = 128;
t        = (1:M.T)/7;
[dYdP,Y] = spm_diff(@(P,M,U)spm_COVID_gen(P,M,U),Ep,M,4,1);

% conditional covariances
%--------------------------------------------------------------------------
Ny    = size(Y,2);
for i = 1:size(Y,2)
    for j = 1:size(dYdP,2)
        D{j} = dYdP{j}(:,i);
    end
    dydp{i}  = spm_cat(D);
    C{i}     = dydp{i}*Cp*dydp{i}';
end

% cumulative death rates
%--------------------------------------------------------------------------
S     = cumsum(Y(:,1));
dSdP  = cumsum(dydp{1});
CS    = dSdP*Cp*dSdP';

% graphics
%--------------------------------------------------------------------------
for i = 1:Ny
    subplot(Ny,2,(i - 1)*2 + 1), hold off
    spm_plot_ci(Y(:,i)',C{i},t), hold on
    try, plot(t(1:numel(Z(:,i))),Z(:,i),'.k'), end
    xlabel('time (weeks)'),ylabel('number of cases/day')
    title(str.outcome{i},'FontSize',16)
    axis square, box off
end

subplot(2,2,2), hold off
spm_plot_ci(S',CS,t), hold on
try, plot(t(1:numel(Z(:,1))),cumsum(Z(:,1)),'.k'), end
xlabel('time (weeks)'),ylabel('number of cases')
title('Cumulative deaths','FontSize',16)
axis square, box off

subplot(2,2,4), hold off
spm_plot_ci(Ep,Cp,[],[],'exp'),               hold on
bar(exp(spm_vec(pE)),1/4,'Edgecolor','none'), hold off
set(gca,'yLim',[0 32])
set(gca,'XTick',1:spm_length(Ep),'Xticklabel',str.names)
ylabel('Parameters','FontSize',16), box off
camorbit(90,0)

