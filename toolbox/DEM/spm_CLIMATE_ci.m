function [Y,C] = spm_CLIMATE_ci(Ep,Cp,Z,U,M,NPI)
% Graphics for climate simulations - with confidence intervals
% FORMAT [Y,C] = spm_CLIMATE_ci(Ep,Cp,Z,U,M,NPI)
% Ep     - posterior expectations
% Cp     - posterior covariances
% Z      - optional empirical data
% U      - indices of outcome
% M      - model
% NPI    - intervention array
%
% Y      - posterior expectation of outcomes
% C      - posterior covariances of outcomes
%
% This routine evaluates a trajectory of outcome variables from a dynamic
% causal model and plots the expected trajectory and accompanying Bayesian
% credible intervals (of 90%). If empirical data are supplied, these will
% be overlaid on the confidence intervals.
%
% Although the DCM is non-linear in the parameters, one can use a
% first-order Taylor expansion to evaluate the confidence intervals in
% terms of how the outcomes change with parameters. This, in combination
% with the well-known overconfidence of variational inference, usually
% requires a slight inflation of uncertainty. Here, the posterior
% covariance is multiplied by a factor of four.
%
% For computational expediency, the confidence intervals are usually
% evaluated as a proportion of the expected value. To evaluate the
% confidence intervals properly, set the global variable CIPLOT to 'true'
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% preliminaries
%--------------------------------------------------------------------------
global CIPLOT, if isempty(CIPLOT); CIPLOT = false; end
if nargin < 6, NPI = []; end

% priors and names for plotting
%--------------------------------------------------------------------------
[~,~,str] = spm_CLIMATE_priors;

% compensate for (variational) overconfidence
%--------------------------------------------------------------------------
Cp = Cp*4;

% evaluate confidence intervals (using a Taylor expansion)
%==========================================================================

% partial derivatives
%----------------------------------------------------------------------
[Y,~,T] = spm_CLIMATE_gen(Ep,M,U,NPI);

% conditional covariances
%----------------------------------------------------------------------
if CIPLOT
    dYdP = spm_diff(@(P,M,U,NPI)spm_CLIMATE_gen(P,M,U,NPI),Ep,M,U,NPI,1);
    C    = dYdP*Cp*dYdP';
else
    C    = diag(abs(Y)/4098);
end

% graphics
%==========================================================================
if numel(Z)
    eplot = min(Z(U).Y) >= 0;
else
    eplot = 1;
end
if eplot
    spm_plot_ci(Y',C,T,[],'exp'), hold on
    if numel(Z)
        plot(Z(U).date,exp(Z(U).Y),'.k')
    end
else
    spm_plot_ci(Y',C,T), hold on
    if numel(Z)
        plot(Z(U).date,Z(U).Y,'.k')
    end
end

datetick('x','yyyy','keeplimits','keepticks')
ylabel('outcome'), xlabel('year')
title(str.outcome(U),'FontSize',14)
box off, spm_axis tight
drawnow
