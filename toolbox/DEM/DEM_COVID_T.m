function [DCM,GCM] = DEM_COVID_T(country,data)
% FORMAT [DCM,GCM] = DEM_COVID(country,data)
% data    - data    to model [default: data = DATA_COVID_JHU]
% country - country to model [default: 'United Kingdom')
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% This routine illustrates the Bayesian model inversion of a generative
% model of coronavirus spread using variational techniques (variational
% Laplace). It illustrates hierarchical Bayesian modelling by first
% inverting a generative model of each country, and then combining the
% posterior densities over the model parameters using parametric empirical
% Bayes to leverage systematic differences between countries, as
% characterised by their population, geographical location etc.
%
% Each subsection produces one or two figures that are described in the
% annotated (Matlab) code. These subsections core various subroutines that
% provide a more detailed description of things like the generative model,
% its priors and the evaluation confidence intervals.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DEM_COVID_T.m 7843 2020-04-30 09:04:45Z karl $

% F: -1.5701e+04 social distancing based upon P(infected)
% F: -1.5969e+04 social distancing based upon P(symptomatic)
% F: -1.5909e+04 social distancing based upon P(waiting)
% F = 0; for i = 1:numel(DCM), F = F + DCM{1}.F; end, disp(F)



% Get data (see DATA_COVID): an array with a structure for each country
%==========================================================================
Y     = DATA_COVID_UK;

% Inversion (i.e., fitting) of empirical data
%==========================================================================

% assemble (Gaussian) priors over model parameters
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% priors for this analysis (use size of population and estimate resistance)
%--------------------------------------------------------------------------
pE.N    = log(66);
pE.r    = log(1/4);

pC.N    = 0;
pC.r    = 1/2;

% optimise parameters for testing (baseline and exponent)
%--------------------------------------------------------------------------
pC.bas  = 1/4;
pC.tes  = 1/4;

% variational Laplace (estimating log evidence (F) and posteriors)
%==========================================================================
% complete model specification
%--------------------------------------------------------------------------
M.G   = @spm_COVID_gen;           % generative function
M.FS  = @(Y)sqrt(Y);              % feature selection  (link function)
M.pE  = pE;                       % prior expectations (parameters)
M.pC  = pC;                       % prior covariances  (parameters)
M.hE  = 0;                        % prior expectation  (log-precision)
M.hC  = 1/64;                     % prior covariances  (log-precision)
M.T   = size(Y,1);                % number of samples
U     = [1 2 6];                  % number of response variables

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);

% assemble prior and posterior estimates (and log evidence)
%--------------------------------------------------------------------------
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.Ep   = Ep;
DCM.Cp   = Cp;
DCM.F    = F;
DCM.Y    = Y;


% save
%--------------------------------------------------------------------------
clear ans
save COVID_UK


% Country specific predictions
%==========================================================================
M.T     = 180;                                    % six-month period
Y       = DCM.Y;                                  % empirical data
Ep      = DCM.Ep;                                 % posterior expectations
Cp      = DCM.Cp;                                 % posterior covariances

% show projections in terms of confidence intervals and superimpose data
%--------------------------------------------------------------------------
spm_figure('GetWin','UK'); clf;
%--------------------------------------------------------------------------
% (predicted outcomes). This figure provides an example of predicted new
% deaths and cases (and recovery rates and its critical care unit
% occupancy) for an exemplar country; here, the United Kingdom. The panels
% on the left shows the predicted outcomes as a function of weeks. The blue
% line corresponds to the expected trajectory, while the shaded areas are
% 90% Bayesian credible intervals. The black dots represent empirical data,
% upon which the parameter estimates are based. The lower right panel shows
% the parameter estimates for the country in question. As in previous
% figures, the prior expectations are shown as in bars over the posterior
% expectations (and credible intervals). The upper right panel illustrates
% the equivalent expectations in terms of cumulative deaths. The key point
% to take from this figure is the quantification of uncertainty inherent in
% the credible or confidence intervals. In other words, uncertainty about
% the parameters propagates through to uncertainty in predicted outcomes.
% This uncertainty changes over time because of the non-linear relationship
% between model parameters and ensemble dynamics. By model design, one can
% be certain about the final states; however, uncertainty about cumulative
% death rates itself accumulates. The mapping from parameters, through
% ensemble dynamics to outcomes is mediated by latent or hidden states. The
% trajectory of these states is illustrated in the subsequent figure.
%--------------------------------------------------------------------------

U   = [1 2 6];
spm_COVID_ci(Ep,Cp,Y,U)

% add seasonal flu rates
%--------------------------------------------------------------------------
FLU = [1692,28330];            % death rate for seasonal flu (per season)
subplot(2,2,2), hold on
x   = get(gca,'XLim');
plot(x,[FLU(1) FLU(1)],'-.r',x,[FLU(2) FLU(2)],'-.r')
spm_axis tight

% and plot latent or hidden states
%--------------------------------------------------------------------------
spm_figure('GetWin','Predictions: UK'); clf;
%--------------------------------------------------------------------------
% (latent causes of observed consequences). The upper panels reproduce the
% expected trajectories of the previous figure, for an example country
% (here the United Kingdom) in. Here, the expected death rate is shown in
% blue, new cases in red, predicted recovery rate in orange and CCU
% occupancy in puce. The black dots correspond to empirical data. The lower
% four panels show the evolution of latent (ensemble) dynamics, in terms of
% the expected probability of being in various states. The first (location)
% panel shows that after about 5 to 6 weeks, there is sufficient evidence
% for the onset of an episode to induce social distancing, such that the
% probability of being found at work falls, over a couple of weeks to
% negligible levels. At this time, the number of infected people increases
% (to about 30%) with a concomitant probability of being infectious a few
% days later. During this time, the probability of becoming immune
% increases monotonically and saturates, within this cell, at about 20
% weeks. Clinically, the probability of becoming symptomatic rises to about
% 20%, with a small probability of developing acute respiratory distress
% and, possibly death. In terms of testing, there is a progressive increase
% in the number of people tested, with a concomitant decrease in those
% untested or waiting for their results. Interestingly, initially the
% number of negative tests increases monotonically, while the proportion of
% positive tests catches up during the peak of the episode. Under these
% parameters, the entire episode lasts for about 12 weeks or three months.
% The increase in (herd) immunity is interesting and will become important
% later. One might ask to what extent these trajectories depend upon
% different model parameters. This is quantified in the next figure.


[Z,X] = spm_COVID_gen(Ep,M,U);
spm_COVID_plot(Z,X,Y,2000,U)

% Sensitivity analysis: which factors determine cumulative deaths?
%==========================================================================
spm_figure('GetWin','Sensitivity: UK'); clf
%--------------------------------------------------------------------------
% (sensitivity analysis). These panels show the change in outcome measures
% (upper panel: death rate. lower panel: new cases). The bar charts are the
% derivatives of outcomes with respect to each of the parameters. Positive
% values (on the right) exacerbate new cases when increased, while,
% conversely, negative values (on the left) decrease new cases. As one
% might expect, increasing social distancing, bed availability and the
% probability of survival outside critical care, tend to decrease death
% rate. Interestingly, increasing both the period of symptoms and ARDS
% decreases overall death rate, because there is more time to recover to an
% asymptomatic state. The lower panel shows the second order derivatives or
% sensitivity. The next figure focuses on the effects of social distancing
% as a way of ameliorating the impact on deaths.

% names{18} = 'testing (post)'; %**
% names{19} = 'testing (initial)'; %**
% names{20} = 'test delay'; %**
% names{21} = 'test selectivity'; %**
% names{22} = 'test exponent'; %**
% names{23} = 'immune period'; %**
% names{24} = 'resistant'; %**

% sensitivity analysis in terms of partial derivatives
%--------------------------------------------------------------------------
M.T = 365*2;
Np  = spm_length(Ep);

p      = [18 19 20 21];
V      = sparse(p,p,1,Np,Np); V = V(:,p);
[dY,P] = spm_diff(@(P,M,U)spm_COVID_gen(P,M,U),Ep,M,1,1,{V});
t      = (1:size(dY,1))/7;
T      = size(Y,1);

subplot(2,1,1)
plot(t,dY,[T,T]/7,[min(dY(:)),max(dY(:))],'-.')
xlabel('time (weeks)')
ylabel('First-order sensitivity'), box off
title('Effect on death rates','FontSize',16),
spm_axis tight, box off
legend(str.names{p}), legend boxoff

% cumulative effects over timea
%--------------------------------------------------------------------------
subplot(2,2,3)
bar(sum(dY(1:T,:)))
xlabel('testing parameter')
ylabel('cumulative deaths')
title('Short-term effect','FontSize',16),
axis square, box off

subplot(2,2,4)
bar(sum(dY(1:end,:)))
xlabel('testing parameter')
ylabel('cumulative deaths')
title('Long-term effect','FontSize',16),
axis square, box off




% Illustrate the effect of social distancing
%==========================================================================
spm_figure('GetWin','Testing: UK'); clf;
%--------------------------------------------------------------------------
% (the effects of social distancing). This figure uses the same format as
% Figure 9. However, here trajectories are reproduced under different
% levels of social distancing; from zero through to 4 (in 16 steps). This
% parameter is the threshold applied to the probability of not being
% infected. In other words, it scores the sensitivity of social distancing
% to the prevalence of the virus in the population. In this example (based
% upon posterior expectations for the United Kingdom), death rates (per
% day) and underlying latent states of the population decrease
% progressively with social distancing. The cumulative death rate is shown
% as a function of social distancing in the upper right panel. The vertical
% line corresponds to the posterior expectation of the social distancing
% threshold for this country. In the next figure, we repeat this analysis
% but looking at the effect of herd immunity.

% P.bas = 1/5000;               % testing sensitivity (to immunity)
% P.sen = 1/10000;              % testing sensitivity (to infection)
% P.del = 2;                    % delay (testing)
% P.tes = 2;                    % test selectivity (for infection)
% P.exp = 1/4;                  % testing exponent


% increase social distancing threshold from 0 to 4
%--------------------------------------------------------------------------
M.T   = 365;                                % one year
P     = Ep;                                 % expansion point
par   = linspace(-4,4,8);                   % range of social distancing
S     = par;
for i = 1:numel(par)
    
    % social distancing threshold
    %----------------------------------------------------------------------
    P.tes = Ep.tes + par(i);
    [Y,X] = spm_COVID_gen(P,M,U);
    S(i)  = sum(Y(:,1));
    
    % plot results and hold graph
    %----------------------------------------------------------------------
    spm_COVID_plot(Y,X,[],[],U)
    for j = 1:6
        subplot(3,2,j), hold on
        set(gca,'ColorOrderIndex',1);
    end
end

% cumulative deaths as a function of social distancing
%--------------------------------------------------------------------------
subplot(3,2,2), hold off
plot(par,S,[1 1]*exp(Ep.sde),[min(S) max(S)],'-.')
title('Social distancing','FontSize',16),
xlabel('social distancing threshold')
ylabel('cumulative deaths')
axis square,box off

disp('lifes saved'), disp(max(S) - min(S))

return
