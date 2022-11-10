function [S,CS,Y,C] = spm_COVID_ci(Ep,Cp,Z,U,M)
% Graphics for coronavirus simulations - with confidence intervals
% FORMAT [S,CS,Y,C] = spm_COVID_ci(Ep,Cp,Z,U,M)
% Ep     - posterior expectations
% Cp     - posterior covariances
% Z      - optional empirical data
% U      - outcomes to evaluate [default: 1:3]
% M      - model
%
% S      - posterior expectation of cumulative outcomes
% CS     - posterior covariances of cumulative outcomes
% Y      - posterior expectation of outcomes
% C      - posterior covariances of outcomes
%
% This routine evaluates a trajectory of outcome variables from a COVID
% model and plots the expected trajectory and accompanying Bayesian
% credible intervals (of 90%). If empirical data are supplied, these will
% be overlaid on the confidence intervals. By default, 128 days
% are evaluated. In addition, posterior and prior expectations are provided
% in a panel.
%
% A single panel is plotted if one output in U is specified
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

% default: number of outcomes to evaluate
%--------------------------------------------------------------------------
if nargin < 4, U = 1:3; end

% priors and names for plotting
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% compensate for (variational) overconfidence
%--------------------------------------------------------------------------
Cp = Cp*4;

% evaluate confidence intervals (using a Taylor expansion)
%==========================================================================

% changes in outcomes with respect to parameters
%--------------------------------------------------------------------------
try, M.T; catch, M.T = 180; end
[dYdP,Y] = spm_diff(@(P,M,U)spm_COVID_gen(P,M,U),Ep,M,U,1);


% conditional covariances
%--------------------------------------------------------------------------
Ny    = size(Y,2);
for i = 1:Ny
    if iscell(dYdP)
        for j = 1:size(dYdP,2)
            D{j} = dYdP{j}(:,i);
        end
        dydp{i}  = spm_cat(D);
    else
        dydp{i}  = dYdP;
    end
    C{i}     = dydp{i}*Cp*dydp{i}';
end

% cumulative death rates
%--------------------------------------------------------------------------
S     = cumsum(Y(:,1));
dSdP  = cumsum(dydp{1});
CS    = dSdP*Cp*dSdP';


% graphics
%==========================================================================
outcome  = str.outcome(U);
if isfield(M,'date')
    t  = (1:M.T) + datenum(M.date,'dd-mm-yyyy');
else
    t  = (1:M.T)/7;
end

% single outcome
%--------------------------------------------------------------------------
if numel(U) == 1
    
    subplot(2,1,1)
    spm_plot_ci(Y(:,i)',C{i},t), hold on
    try, plot(t(1:numel(Z(:,i))),Z(:,i),'.k'), end
    ylabel('number of cases/day')
    title(outcome,'FontSize',16)
    
    
    % label time
    %----------------------------------------------------------------------
    if isfield(M,'date')
        datetick('x','mmmdd')
        xlabel('date')
    else
        xlabel('time (weeks)')
    end
    
    box off, spm_axis tight
    YLim = get(gca,'YLim'); YLim(1) = 0; set(gca,'YLim',YLim);
    return
end

% graphics
%--------------------------------------------------------------------------
for i = 1:Ny
    subplot(Ny,2,(i - 1)*2 + 1)
    spm_plot_ci(Y(:,i)',C{i},t), hold on
    try, plot(t(1:numel(Z(:,i))),Z(:,i),'.k'), end
    ylabel('number of cases/day')
    title(outcome{i},'FontSize',16)
    
    % label time
    %----------------------------------------------------------------------
    if isfield(M,'date')
        datetick('x','mmmdd')
        xlabel('date')
    else
        xlabel('time (weeks)')
    end
    axis square, box off, spm_axis tight
    YLim = get(gca,'YLim'); YLim(1) = 0; set(gca,'YLim',YLim);
end

subplot(2,2,2)
spm_plot_ci(S',CS,t), hold on
try, plot(t(1:numel(Z(:,1))),cumsum(Z(:,1)),'.k'), end
xlabel('time (weeks)'),ylabel('number of cases')
title('Cumulative numbers','FontSize',16)

% label time
%----------------------------------------------------------------------
if isfield(M,'date')
    datetick('x','mmmdd')
    xlabel('date')
else
    xlabel('time (weeks)')
end
axis square, box off, spm_axis tight
YLim = get(gca,'YLim'); YLim(1) = 0; set(gca,'YLim',YLim);


subplot(2,2,4)
spm_plot_ci(Ep,Cp,[],[],'exp'),               hold on
bar(exp(spm_vec(pE)),1/4,'Edgecolor','none'), hold off
set(gca,'XTick',1:spm_length(Ep),'Xticklabel',str.names)
ylabel('Parameters','FontSize',16), box off
camorbit(90,0)

