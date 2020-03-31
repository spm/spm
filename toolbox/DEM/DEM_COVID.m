function [DCM,GCM] = DEM_COVID(data)
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% This routine illustrates the Bayesian model inversion of a generative
% model of coronavirus spread using variational techniques (variational
% Laplace). It illustrates hierarchical Bayesian modelling by first
% inverting a generative model of each country, and then combining the
% posterior densities over the model parameters using parametric empirical
% Bayesto leverage systematic differences between countries, as
% characterised by their population, geographical location etc.
%
% This routine produces a series of figures illustrating parameter
% estimates and their associated confidence intervals. It then moves on to
% look at timeseries predictions for a specific country (currently
% hardcoded to be United Kingdom). It considers the effects of various
% model parameters such as those encoding the degree of social distancing
% and initial (herd) immunity.
%
% Each subsection produces one or two figures that are described in the
% annotated (Matlab)code. These subsections core various subroutines that
% provide a more detailed description of things like the generative model,
% its priors and the evaluation confidence intervals.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_COVID.m 7809 2020-03-31 11:55:09Z karl $

% F: -1.5701e+04 social distancing based upon P(infected)
% F: -1.5969e+04 social distancing based upon P(symptomatic)
% F: -1.5909e+04 social distancing based upon P(waiting)
% F = 0; for i = 1:numel(DCM), F = F + DCM{1}.F; end, disp(F)



% Get data (see DATA_COVID): an array with a structure for each country
%==========================================================================
if nargin < 1,  data = DATA_COVID_JHU; end


% Inversion (i.e., fitting) of empirical data
%==========================================================================

% assemble (Gaussian) priors over model parameters
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% Bayesian inversion (placing posteriors in a cell array of structures)
%--------------------------------------------------------------------------
GCM   = cell(size(data(:)));
for i = 1:numel(data)
    
    % data for this country (here, and positive test rates)
    %----------------------------------------------------------------------
    set(gcf,'name',data(i).country)
    Y   = [data(i).death, data(i).cases];
    
    % variational Laplace (estimating log evidence (F) and posteriors)
    %======================================================================
    [F,Ep,Cp,pE,pC] = spm_COVID(Y,pE,pC);
    
    % assemble prior and posterior estimates (and log evidence)
    %----------------------------------------------------------------------
    DCM.M.pE = pE;
    DCM.M.pC = pC;
    DCM.Ep   = Ep;
    DCM.Cp   = Cp;
    DCM.F    = F;
    
    % save this country (in a cell array)
    %----------------------------------------------------------------------
    GCM{i}   = DCM;
    
end

% Between country analysis (hierarchical or parametric empirical Bayes)
%==========================================================================
spm_figure('GetWin','Bayesian model reduction'); clf;
%--------------------------------------------------------------------------
% (Bayesian model comparison). This figure with reports the result of
% Bayesian model comparison (a.k.a. Bayesian model reduction). In this
% instance, the models compared are at the second or between country level.
% In other words, the models compared contained all combinations of (second
% level) parameters (a parameter is removed by setting its prior covariance
% to zero). If the model evidence increases – in virtue of reducing model
% complexity – then this parameter is redundant. The redundant parameters
% are shown in the lower panels by comparing the posterior expectations
% before and after Bayesian model reduction. The blue bars correspond to
% posterior expectations, while the pink bars denote 90% Bayesian credible
% intervals. The key thing to take from this analysis is that a large
% number of second level parameters are redundant. These second level
% parameters encode the effects of population size and geographical
% location, on each of the parameters of the generative model. The next
% figure illustrates the nonredundant effects that can be inferred with
% almost a 100% posterior confidence.

% general linear model (constant, log-population, latitude and longitude)
%--------------------------------------------------------------------------
lat    = spm_vec([data.lat]);
lon    = spm_vec([data.long]);
lat    = lat*2*pi/360;
lon    = lon*2*pi/360;
X      = [];
Xn     = {'const','log cell size'};
for  i = 1:4
    X  = [X sin(i*lon) sin(i*lat)];
    Xn = {Xn{:}, sprintf('lat(%d)',i), sprintf('lon(%d)',i)};
end

X      = [log(spm_vec([data.pop])) X];
X      = [ones(numel(data),1) X];
X      = spm_orth(X,'norm');
X(:,1) = 1;
GLM.X  = X;
GLM.Xnames = Xn;

% parametric empirical Bayes (with random effects in str.field)
%--------------------------------------------------------------------------
[PEB,DCM] = spm_dcm_peb(GCM,GLM,str.field);

% Bayesian model averaging (over reduced models), testing for GLM effects
%--------------------------------------------------------------------------
[BMA,BMR] = spm_dcm_bmr_all(PEB,str.field);

% Bayesian parameter averaging (over countries)
%--------------------------------------------------------------------------
BPA       = spm_dcm_bpa(DCM,'nocd');


% illustrate the largest between country effects
%==========================================================================
spm_figure('GetWin','Second level effects'); clf;
%--------------------------------------------------------------------------
% (between country effects). This figure shows the relationship between
% certain parameters of the generative model and the explanatory variables
% in a general linear model of between country effects. The examples are
% based upon a ranking of the absolute value of the second level parameter;
% namely, the contribution of an explanatory variable to a model parameter.
% Here, the effective size of the population cell appears to depend upon
% the latitude of a country.

% assemble parameters
%--------------------------------------------------------------------------
P     = [];                                   % posterior expectations
C     = [];                                   % posterior variances
for i = 1:numel(DCM)
    P(:,i) = spm_vec(DCM{i}.Ep);
end
P     = P(PEB.Pind,:);
Pname = str.names(PEB.Pind);

% find largest absolute (second level) effects and plot
%--------------------------------------------------------------------------
Ep      = abs(PEB.Ep);
Ep(:,1) = 0;
Sp      = sort(Ep(:),'descend');

for i = 1:6
    [I,J] = find(Ep == Sp(i));
    subplot(3,2,i), plot(X(:,J),P(I,:),'.','MarkerSize',32,'Color',[0.8 0.8 1])
    xlabel(PEB.Xnames{J}),ylabel([ 'log ' Pname{I}])
    title(Pname{I},'FontSize',16),axis square, box off
end

% report Bayesian parameter averages, in relation to priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Parameter estimates'); clf;
%--------------------------------------------------------------------------
% (Bayesian parameter averages). This figure reports the Bayesian parameter
% averages over countries following a hierarchical or parametric empirical
% Bayesian analysis that tests for – and applies shrinkage priors to –
% posterior parameter estimates for each country. The upper panel shows the
% parameters as estimated in log space, while the lower panel shows the
% same results for the corresponding scale parameters (scale parameters are
% nonnegative parameters). The blue bars report posterior expectations,
% while the thin red bars are prior expectations. The pink bars denote 90%
% Bayesian confidence or credible intervals. One can interpret these
% parameters as the average value for any given parameter of the generative
% model, to which a random (country specific) effect is added to generate
% the ensemble dynamics for each country. In turn, these ensemble
% distributions determine the likelihood of various outcome measures under
% larger number (i.e., Poisson) assumptions. For example, the average
% number of interpersonal contacts (that are potentially contagious) at
% home is about four, while at work it is about 24. The two phases of
% infection (infected and contagious) are roughly the same at around five
% days. Notice that these estimates are not based upon empirical data.
% These are the periods that provide the best explanation for the new cases
% and deaths, reported from each country.

Ep = BPA.Ep;                                  % posterior expectations                         
Cp = BPA.Cp;                                  % posterior covariances

subplot(2,1,1)
spm_plot_ci(Ep,Cp), hold on, bar(spm_vec(pE),1/4), hold off
ylabel('log parameters','FontSize',16)
set(gca,'XTick',1:spm_length(Ep),'Xticklabel',str.names)
camorbit(90,0), axis square, box off

subplot(2,1,2)
spm_plot_ci(Ep,Cp,[],[],'exp')
set(gca,'yLim',[0 32])
set(gca,'XTick',1:spm_length(Ep),'Xticklabel',str.names)
ylabel('Parameters','FontSize',16)
camorbit(90,0), axis square, box off


% Differences among countries, in terms of parameters
%==========================================================================
spm_figure('GetWin','Parameters');
%--------------------------------------------------------------------------
% (differences among countries). This figure reports the differences among
% countries in terms of selected parameters of the generative model,
% ranging from the size of a cell (i.e., the effective size of an infected
% population), through to the probability of dying when in critical care.
% Interesting country specific differences here include an apparent
% attenuation of social distancing responses, relative to other countries,
% in the United States and Australia. The blue bars represent the posterior
% expectations, while the pink bars are 90% Bayesian credible intervals.
% Notice that these intervals are not symmetrical about the mean because
% scale parameters are plotted here – as opposed to the log parameters. The
% next figure illustrates the predictions – in terms of new deaths and
% cases – based upon these parameter estimates.

% assemble parameters
%--------------------------------------------------------------------------
P     = [];                                   % posterior expectations
C     = [];                                   % posterior variances
for i = 1:numel(DCM)
    P(:,i) = spm_vec(DCM{i}.Ep);
    C(:,i) = diag(DCM{i}.Cp);
end

% report selected parameters (see spm_COVID_priors)
%--------------------------------------------------------------------------
Np    = spm_length(pE);
j     = 1:Np; j([1,3,5,7,11,15,18,19,21]) = [];
for i = 1:length(j)
    subplot(4,3,i) 
    spm_plot_ci(P(j(i),:)',C(j(i),:),[],[],'exp')
    title(str.names{j(i)},'FontSize',16)
    xlabel('country'), axis square, box off
end


% Country specific predictions
%==========================================================================
country = 'United Kingdom';                       % country to predict
i       = find(ismember({data.country},country)); % country index
Y       = [data(i).death, data(i).cases];         % empirical data
M.T     = 180;                                    % six-month period
Ep      = DCM{i}.Ep;                              % posterior expectations
Cp      = DCM{i}.Cp;                              % posterior covariances

% show projections in terms of confidence intervals and superimpose data
%--------------------------------------------------------------------------
spm_figure('GetWin',country); clf;
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

spm_COVID_ci(Ep,Cp,Y)

% and plot latent or hidden states
%--------------------------------------------------------------------------
spm_figure('GetWin',['Predictions: ' country]); clf;
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
% (to about 18%) with a concomitant probability of being infectious a few
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

[Z,X] = spm_COVID_gen(DCM{i}.Ep,M,4);
spm_COVID_plot(Z,X,Y)


% Sensitivity analysis: which factors determine cumulative deaths?
%==========================================================================
spm_figure('GetWin',['Sensitivity: ' country]);
%--------------------------------------------------------------------------
% (sensitivity analysis). These panels show the change in outcome measures
% (upper panel: death rate. lower panel: new cases). The bar charts are the
% derivatives of outcomes with respect to each of the parameters. Positive
% values (on the right) exacerbate new cases when increased, while,
% conversely, negative values (on the left) decrease new cases. As one
% might expect, increasing social distancing, bed availability and the
% probability of survival outside critical care, tend to decrease death
% rate. Interestingly, increasing both the period of symptoms and ARDS
% decreases overall death rate; presumably, because there is more time to
% recover to an asymptomatic state and the probability of infecting someone
% else – when symptomatic – is attenuated. The next figure focuses on the
% effects of social distancing as a way of ameliorating the impact on
% deaths.

% sensitivity analysis in terms of partial derivatives
%--------------------------------------------------------------------------
[dYdP,Y] = spm_diff(@(P,M,U)spm_COVID_gen(P,M,U),Ep,M,2,1);
Ny       = size(Y,2);

% cumulative effects over time
%--------------------------------------------------------------------------
for i = 1:Ny
    for j = 1:size(dYdP,2)
        D{j} = dYdP{j}(:,i);
    end
    dRdP{i}  = sum(spm_cat(D));
end

% plot results
%--------------------------------------------------------------------------
for i = 1:Ny
    subplot(Ny,1,i)
    bar(dRdP{i})
    set(gca,'XTick',1:spm_length(Ep),'Xticklabel',str.names,'FontSize',8)
    ylabel(str.outcome{i},'FontSize',16), box off
    camorbit(90,0),axis square
end


% Illustrate the effect of social distancing
%==========================================================================
spm_figure('GetWin',['Social distancing:' country]); clf;
%--------------------------------------------------------------------------
% (the effects of social distancing). This figure uses the same format as
% Figure 9. However, here trajectories are reproduced under different
% levels of social distancing; from zero through to 4 (in 16 steps). This
% parameter is the exponent applied to the probability of not being
% infected. In other words, it scores the sensitivity of social distancing
% to the prevalence of the virus in the population. In this example (based
% upon posterior expectations for the United Kingdom), death rates (per
% day) and underlying latent states of the population decrease
% progressively with social distancing. The cumulative death rate is shown
% as a function of social distancing in the upper right panel. The vertical
% line corresponds to the posterior expectation of the social distancing
% exponent for this country. These results suggest that social distancing
% relieves pressure on critical care capacities and ameliorates cumulative
% deaths by about 3700 people. This is roughly four times the number of
% people who die in the equivalent period due to road traffic accidents. In
% the next figure, we repeat this analysis but looking at the effect of
% herd immunity.

% increase social distancing exponent from 0 to 4
%--------------------------------------------------------------------------
P     = Ep;                                    % expansion point
sde   = linspace(0,4,16);                      % range of social distancing
S     = sde;
for i = 1:numel(sde)
    P.sde  = Ep.sde + log(sde(i) + 1e-6);
    [Y,X]  = spm_COVID_gen(P,M,1);
    S(i)   = sum(Y);
    
    % plot results and hold graph
    %----------------------------------------------------------------------
    spm_COVID_plot(Y,X)
    for j = 1:6, subplot(3,2,j), hold on, end
end

% cumulative deaths as a function of social distancing
%--------------------------------------------------------------------------
subplot(3,2,2), hold off
plot(sde,S,[1 1]*exp(Ep.sde),[min(S) max(S)],'-.')
title('Social distancing','FontSize',16), 
xlabel('social distancing exponent')
ylabel('cumulative deaths')
axis square,box off


% Herd immunity
%==========================================================================
spm_figure('GetWin',['Herd immunity:' country]); clf;
%--------------------------------------------------------------------------
% (herd immunity). This figure reproduces the format of the previous
% figure. However, here, we increased the initial proportion of the cell
% (i.e., population) who were initially immune. Increasing the initial
% immunity dramatically decreases death rates with a fall in the cumulative
% deaths from several thousand to negligible levels with a herd immunity of
% about 70%. The dashed line in the upper panel shows the equivalent deaths
% over the same time period due to seasonal flu (based upon 2019 figures).
% This death rate would require an initial or herd immunity of about 60%.
% It is interesting to return to Figure 6 and identify at what point –
% during the course of the infection episode – this level of herd immunity
% is obtained.

%--------------------------------------------------------------------------
% Public Health England estimates that on average 17,000 people have died
% from the flu in England annually between 2014/15 and 2018/19. However,
% the yearly deaths vary widely, from a high of 28,330 in 2014/15 to a low
% of 1,692 in 2018/19. Public Health England does not publish a mortality
% rate for the flu.

% The Department for Transport (DfT) has announced there were 1,784
% reported road deaths in 2018, compared to 1,793 reported in 2017 – a 1%
% fall. There were 25,511 people seriously injured in reported road traffic
% accidents in 2018, compared to 24,831 in 2017 – a 3% year-on-year
% increase
%--------------------------------------------------------------------------

% progressively increase initial immunity
%--------------------------------------------------------------------------
FLU   = 1692/365;                  % death rate for seasonal flu (per day)
m     = linspace(0,1,16);          % range of initial immunity
P     = Ep;
S     = m;
for i = 1:numel(m)
    P.m   = log(m(i) + 1e-6);
    [Y,X] = spm_COVID_gen(P,M,1);
    S(i)  = sum(Y);
    
    % plot results and hold graph
    %----------------------------------------------------------------------
    spm_COVID_plot(Y,X)
    for j = 1:6, subplot(3,2,j), hold on, end
    
end

% plot
%--------------------------------------------------------------------------
subplot(3,2,2), hold off
plot(m,S,m,FLU*M.T*m.^0,'-.')
title('Herd immunity','FontSize',16), 
xlabel('proportion immune')
ylabel('cumulative deaths')
axis square,box off


% demonstrate routines: predictive validity
%==========================================================================
% (predictive validity – early). This figure uses the same format as Figure
% 5; however, here, the posterior estimates are based upon partial data,
% from early in the timeseries for an exemplar country. These estimates are
% based upon the empirical priors following parametric empirical Bayes. The
% red dots show the outcomes that were observed but not used to estimate
% the expected trajectories (or confidence intervals). This example
% illustrates the predictive validity of the estimates for a two-week
% period following the last datapoint. This captures the rise to the peak
% of new cases in Italy.

% remove ( > T) data from country ( = i)
%--------------------------------------------------------------------------
% T           = 35;                         % number of days to withhold
% i           = 1;                          % country index
T             = 14;                         % number of days to withhold
i             = 2;                          % country index



% use priors from parametric empirical Bayes
%--------------------------------------------------------------------------
pE            = DCM{i}.M.pE;
pC            = DCM{i}.M.pC;
data          = DATA_COVID_JHU;
data(i).death = data(i).death(1:end - T);
data(i).cases = data(i).cases(1:end - T);

% invert (using incomplete data) and plot confidence intervals
%--------------------------------------------------------------------------
Y             = [data(i).death, data(i).cases];
[F,Ep,Cp]     = spm_COVID(Y,pE,pC);

spm_figure('GetWin','predictive validity'); clf
spm_COVID_ci(Ep,Cp,Y)

% retrieve and overlay withheld data
%--------------------------------------------------------------------------
data = DATA_COVID_JHU;
Y    = [data(i).death, data(i).cases];
NY   = size(Y,1);
t    = (1:NY)/7;
CY   = cumsum(Y(:,1));
i    = (NY - T):NY;
t    = t(i);

spm_figure('GetWin','predictive validity');
subplot(4,2,1), hold on, plot(t,Y(i,1),'.r','MarkerSize',16)
subplot(4,2,3), hold on, plot(t,Y(i,2),'.r','MarkerSize',16)
subplot(2,2,2), hold on, plot(t,CY(i), '.r','MarkerSize',16)

return



% auxiliary routines
%__________________________________________________________________________


% demonstrate routines: face validation of inversion scheme
%==========================================================================

% Gaussian priors over model parameters
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% Sample model parameters from prior distribution
%--------------------------------------------------------------------------
dP    = diag(sqrt(spm_vec(pC)))*randn(spm_length(pE),1);
P     = spm_unvec(spm_vec(pE) + dP,pE);
[Y,X] = spm_COVID_gen(P,M,2);

% plot synthetic data
%--------------------------------------------------------------------------
spm_COVID_plot(Y,X)

% Variational Laplace
%==========================================================================
[F,Ep,Cp,pE,pC,Eh] = spm_COVID(Y,pE,pC);

% compare true and estimated model parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Parameter estimates'); clf;

subplot(2,1,1)
spm_plot_ci(Ep,Cp), hold on, bar(spm_vec(pE),1/4), hold off
ylabel('log parameters','FontSize',16)
set(gca,'XTick',1:spm_length(Ep),'Xticklabel',str.names)
camorbit(90,0), axis square

subplot(2,1,2)
spm_plot_ci(Ep,Cp,[],[],'exp')
set(gca,'yLim',[0 32])
set(gca,'XTick',1:spm_length(Ep),'Xticklabel',str.names)
ylabel('Parameters','FontSize',16)
camorbit(90,0), axis square

% plot data fit
%--------------------------------------------------------------------------
spm_figure('GetWin','Predictions');
[Z,X] = spm_COVID_gen(Ep,M,4);
spm_COVID_plot(Z,X,Y)


