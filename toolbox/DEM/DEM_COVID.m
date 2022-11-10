function [DCM,GCM] = DEM_COVID(country,data)
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
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% F: -1.5701e+04 social distancing based upon P(infected)
% F: -1.5969e+04 social distancing based upon P(symptomatic)
% F: -1.5909e+04 social distancing based upon P(waiting)
% F = 0; for i = 1:numel(DCM), F = F + DCM{1}.F; end, disp(F)



% Get data (see DATA_COVID): an array with a structure for each country
%==========================================================================
if nargin < 2, data    = DATA_COVID_JHU(16); end
if nargin < 1, country = 'United Kingdom';   end

% Inversion (i.e., fitting) of empirical data
%==========================================================================
Fsi = spm_figure('GetWin','SI'); clf;

% assemble (Gaussian) priors over model parameters
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;
hC          = 1/64;

% Bayesian inversion (placing posteriors in a cell array of structures)
%--------------------------------------------------------------------------
GCM   = cell(size(data(:)));
for i = 1:numel(data)
    
    % data for this country (here, and positive test rates)
    %----------------------------------------------------------------------
    set(Fsi,'name',data(i).country)
    Y = [data(i).death, data(i).cases];
   
    % variational Laplace (estimating log evidence (F) and posteriors)
    %======================================================================
    [F,Ep,Cp,pE,pC] = spm_COVID(Y,pE,pC,hC);
    
    
    % assemble prior and posterior estimates (and log evidence)
    %----------------------------------------------------------------------
    DCM.M.pE = pE;
    DCM.M.pC = pC;
    DCM.Ep   = Ep;
    DCM.Cp   = Cp;
    DCM.F    = F;
    DCM.Y    = Y;
    
    % save this country (in a cell array)
    %----------------------------------------------------------------------
    GCM{i}   = DCM;
    
end

% Between country analysis (hierarchical or parametric empirical Bayes)
%==========================================================================
spm_figure('GetWin','BMR - all'); clf;
%--------------------------------------------------------------------------
% (Bayesian model comparison). This figure with reports the result of
% Bayesian model comparison (a.k.a. Bayesian model reduction). In this
% instance, the models compared are at the second or between country level.
% In other words, the models compared contained all combinations of (second
% level) parameters (a parameter is removed by setting its prior covariance
% to zero). If the model evidence increases - in virtue of reducing model
% complexity - then this parameter is redundant. The redundant parameters
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
% M.X      - 2nd-level design matrix: X(:,1) = ones(N,1) [default]
% M.bE     - 3rd-level prior expectation [default: DCM{1}.M.pE]
% M.bC     - 3rd-level prior covariance  [default: DCM{1}.M.pC/M.alpha]
% M.pC     - 2nd-level prior covariance  [default: DCM{1}.M.pC/M.beta]
%
% M.alpha  - optional scaling to specify M.bC [default = 1]
% M.beta   - optional scaling to specify M.pC [default = 16]
%--------------------------------------------------------------------------
lat    = spm_vec([data.lat]);
lon    = spm_vec([data.long]);
lat    = lat*2*pi/360;
lon    = lon*2*pi/360;
X      = [];
Xn     = {'const','log(N)'};
for  i = 1:4
    X  = [X sin(i*lon) sin(i*lat)];
    Xn = [Xn(:)', {sprintf('lat(%d)',i)}, {sprintf('lon(%d)',i)}];
end

% design matrix of explanatory variables
%--------------------------------------------------------------------------
X      = [log(spm_vec([data.pop])) X];
X      = [ones(numel(data),1) X];
X      = spm_orth(X,'norm');
X(:,1) = 1;

% place in general linear model
%--------------------------------------------------------------------------
GLM.X      = X;
GLM.Xnames = Xn;

% parametric empirical Bayes (with random effects in str.field)
%==========================================================================
[PEB,DCM] = spm_dcm_peb(GCM,GLM,str.field);

% Bayesian model averaging (over reduced models), testing for GLM effects
%--------------------------------------------------------------------------
[BMA,BMR] = spm_dcm_bmr_all(PEB,str.field);

% Repeat inversion using parametric empirical priors
%==========================================================================
for i = 1:numel(DCM)

    % variational Laplace
    %----------------------------------------------------------------------
    set(Fsi,'name',data(i).country)
    [F,Ep,Cp] = spm_COVID(DCM{i}.Y,DCM{i}.M.pE,DCM{i}.M.pC,hC);
    
    % assemble prior and posterior estimates (and log evidence)
    %----------------------------------------------------------------------
    DCM{i}.Ep = Ep;
    DCM{i}.Cp = Cp;
    DCM{i}.F  = F;
    
end

% Bayesian parameter averaging (over countries)
%--------------------------------------------------------------------------
BPA       = spm_dcm_bpa(DCM,'nocd');

% save
%--------------------------------------------------------------------------
clear Fsi ans
save COVID_DCM


% Illustrate the largest between country effects
%==========================================================================
spm_figure('GetWin','Second level effects'); clf;
%--------------------------------------------------------------------------
% (between country effects). This figure shows the relationship between
% certain parameters of the generative model and the explanatory variables
% in a general linear model of between country effects. The examples are
% based upon a ranking of the absolute value of the second level parameter;
% namely, the contribution of an explanatory variable to a model parameter.
% The lower panel shows the (absolute) parameters in image format

% assemble parameters
%--------------------------------------------------------------------------
P     = [];                                   % posterior expectations
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

for i = 1:4
    [I,J] = find(Ep == Sp(i));
    subplot(3,2,i), plot(X(:,J),P(I,:),'.','MarkerSize',32,'Color',[0.8 0.8 1])
    xlabel(PEB.Xnames{J}),ylabel([ 'log ' Pname{I}])
    title(Pname{I},'FontSize',16),axis square, box off
end

% GLM (second level) parameters
%--------------------------------------------------------------------------
subplot(3,1,3)
i = PEB.Pind;
imagesc(Ep), title('Parameters of GLM','FontSize',16)
set(gca,'XTick',1:numel(GLM.Xnames) ,'Xticklabel',GLM.Xnames)
set(gca,'YTick',1:numel(str.names(i)),'Yticklabel',str.names(i))
try, set(gca,'XTickLabelRotation',90), end
axis square, box off


% report Bayesian parameter averages, in relation to priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Parameter estimates'); clf;
%--------------------------------------------------------------------------
% (Bayesian parameter averages). This figure reports the Bayesian parameter
% averages over countries following a hierarchical or parametric empirical
% Bayesian analysis that tests for - and applies shrinkage priors to -
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
% larger number (i.e., Poisson) assumptions.

Ep = BPA.Ep;                                  % posterior expectations                         
Cp = BPA.Cp;                                  % posterior covariances

subplot(2,1,1)
spm_plot_ci(Ep,Cp), hold on, bar(spm_vec(pE),1/4), hold off
ylabel('log parameters','FontSize',16)
set(gca,'XTick',1:spm_length(Ep),'Xticklabel',str.names)
camorbit(90,0), axis square, box off

subplot(2,1,2)
spm_plot_ci(Ep,Cp,[],[],'exp')
set(gca,'XTick',1:spm_length(Ep),'Xticklabel',str.names)
ylabel('Parameters','FontSize',16)
camorbit(90,0), axis square, box off


% Differences among countries, in terms of parameters
%==========================================================================
spm_figure('GetWin','Parameters'); clf;
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
% scale parameters are plotted here - as opposed to the log parameters. The
% next figure illustrates the predictions - in terms of new deaths and
% cases - based upon these parameter estimates.

% assemble parameters
%--------------------------------------------------------------------------
P     = [];                                   % posterior expectations
C     = [];                                   % posterior variances
for i = 1:numel(DCM)
    P(:,i) = spm_vec(DCM{i}.Ep);
    C(:,i) = diag(DCM{i}.Cp);
end

% names{1}  = 'initial cases'; %**
% names{2}  = 'population size';
% names{3}  = 'initial immunity';
% names{4}  = 'P(work | home)';
% names{5}  = 'social distancing';
% names{6}  = 'bed availability';  
% names{7}  = 'contacts: home';
% names{8}  = 'contacts: work';
% names{9}  = 'transmission strength';
% names{10} = 'infected period';
% names{11} = 'contagious period';
% names{12} = 'incubation period';
% names{13} = 'P(ARDS | symptoms)';
% names{14} = 'symptomatic period';
% names{15} = 'time in CCU';
% names{16} = 'P(fatality | CCU)';
% names{17} = 'P(survival | home)';
% names{18} = 'trace and test'; %**
% names{19} = 'testing latency'; %**
% names{20} = 'test delay'; %**
% names{21} = 'test selectivity'; %**
% names{22} = 'sustained testing'; %**
% names{23} = 'baseline testing'; %**
% names{24} = 'immune period'; %**
% names{25} = 'resistance'; %**

% report selected parameters (see spm_COVID_priors)
%--------------------------------------------------------------------------
p     = 1:size(P,1); p([1 3 17 18 20 21 22 23 24]) = [];
p     = p(1:16);
for i = 1:length(p)
    
    % posterior density
    %----------------------------------------------------------------------
    subplot(4,4,i)
    Ep   = P(p(i),:);
    Cp   = C(p(i),:);
    spm_plot_ci(Ep',Cp,[],[],'exp'), hold on
    title(str.names{p(i)},'FontSize',12)
    xlabel('country'), axis square, box off
    
    % country with greatest map estimate (and United Kingdom)
    %----------------------------------------------------------------------
    [d,j] = max(exp(Ep));
    text(j,d,data(j).country,'FontSize',8);
    [d,j] = min(exp(Ep));
    text(j,d,data(j).country,'FontSize',8);
    j     = find(ismember({data.country},'United Kingdom'));
    text(j,exp(Ep(j)),'*','FontSize',12,'Color','r','HorizontalAlignment','center');
    hold off
    
end


% Country specific predictions
%==========================================================================
i   = find(ismember({data.country},country)); % country index
M   = DCM{i}.M;                               % model
Y   = DCM{i}.Y;                               % empirical data
Ep  = DCM{i}.Ep;                              % posterior expectations
Cp  = DCM{i}.Cp;                              % posterior covariances
M.T = 180;                                    % six-month period

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
%--------------------------------------------------------------------------
% Public Health England estimates that on average 17,000 people have died
% from the flu in England annually between 2014/15 and 2018/19. However,
% the yearly deaths vary widely, from a high of 28,330 in 2014/15 to a low
% of 1,692 in 2018/19. Public Health England does not publish a mortality
% rate for the flu.
%--------------------------------------------------------------------------
spm_COVID_ci(Ep,Cp,Y);

% add seasonal flu rates
%--------------------------------------------------------------------------
FLU = [1692,28330];            % death rate for seasonal flu (per season)
subplot(2,2,2), hold on
x   = get(gca,'XLim');
plot(x,[FLU(1) FLU(1)],'-.r',x,[FLU(2) FLU(2)],'-.r')
spm_axis tight

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
%--------------------------------------------------------------------------
[Z,X] = spm_COVID_gen(Ep,M,1:3);
spm_COVID_plot(Z,X,Y);


% Sensitivity analysis: which factors determine cumulative deaths?
%==========================================================================
spm_figure('GetWin',['Sensitivity: ' country]); clf
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

% sensitivity analysis in terms of partial derivatives
%--------------------------------------------------------------------------
[ddY,dY] = spm_diff(@(P,M,U)spm_COVID_gen(P,M,U),Ep,M,1,[1,1]);
Np       = spm_length(Ep);

% cumulative effects over time
%--------------------------------------------------------------------------
DY    = sum(dY);
for i = 1:Np
    DDY{i} = sum(ddY{i});
end
DDY   = spm_cat(DDY');

% plot results
%--------------------------------------------------------------------------
subplot(2,1,1)
bar(DY)
set(gca,'XTick',1:Np,'Xticklabel',str.names,'FontSize',8)
ylabel('First-order sensitivity','FontSize',16), box off
camorbit(90,0),axis square

subplot(2,1,2)
imagesc(DDY)
set(gca,'YTick',1:Np,'Yticklabel',str.names,'FontSize',8)
title('Second-order sensitivity','FontSize',16), box off
axis square


% Illustrate the effect of social distancing
%==========================================================================
spm_figure('GetWin',['Social distancing:' country]); clf;
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

% increase social distancing threshold from -4 to 4 (log scaling)
%--------------------------------------------------------------------------
P     = Ep;                                 % expansion point
sde   = linspace(-4,4,16);                  % range of social distancing
S     = sde;
for i = 1:numel(sde)
    
    % social distancing threshold
    %----------------------------------------------------------------------
    P.sde = Ep.sde + sde(i);
    [Y,X] = spm_COVID_gen(P,M,1);
    S(i)  = sum(Y(:,1));
    
    % plot results and hold graph
    %----------------------------------------------------------------------
    spm_COVID_plot(Y,X)
    for j = 1:6
        subplot(3,2,j), hold on
        set(gca,'ColorOrderIndex',1);
    end
end

% cumulative deaths as a function of social distancing
%--------------------------------------------------------------------------
subplot(3,2,2), hold off
plot(sde,S,[1 1]*exp(Ep.sde),[min(S) max(S)],'-.')
title('Social distancing','FontSize',16),
xlabel('social distancing threshold')
ylabel('cumulative deaths')
axis square,box off

disp('lifes saved'), disp(max(S) - min(S))
    

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
% over the same time period due to seasonal flu (based upon 2014/19
% figures). This death rate would require an initial or herd immunity of
% about 60%. It is interesting to return to Figure 6 and identify at what
% point - during the course of the infection episode - this level of herd
% immunity is obtained.

%--------------------------------------------------------------------------
% The Department for Transport (DfT) has announced there were 1,784
% reported road deaths in 2018, compared to 1,793 reported in 2017 - a 1%
% fall. There were 25,511 people seriously injured in reported road traffic
% accidents in 2018, compared to 24,831 in 2017 - a 3% year-on-year
% increase
%--------------------------------------------------------------------------

% progressively increase initial immunity
%--------------------------------------------------------------------------
m     = linspace(0,1,16);
P     = Ep;
S     = m;
for i = 1:numel(m)
    P.m   = log(m(i) + 1e-6);
    [Y,X] = spm_COVID_gen(P,M,1);
    S(i)  = sum(Y);
    
    % plot results and hold graph
    %----------------------------------------------------------------------
    spm_COVID_plot(Y,X)
    for j = 1:6
        subplot(3,2,j), hold on
        try, set(gca,'ColorOrderIndex',1); end
    end
    
end

% plot
%--------------------------------------------------------------------------
subplot(3,2,2), hold off
S(S < 0) = 0;
plot(m,S,m,FLU(1)*m.^0,'-.r',m,FLU(2)*m.^0,'-.r')
title('Herd immunity','FontSize',16), 
xlabel('proportion immune')
ylabel('cumulative deaths')
axis square,box off

% ensemble dynamics in terms of basic reproduction rate
%==========================================================================
spm_figure('GetWin',['Reproduction rate:' country]); clf;
%--------------------------------------------------------------------------
% (basic reproduction ratio). This figure plots the predicted death rates
% for the country in question above the concomitant fluctuations
% in the basic reproduction rate (R0) and herd immunity. The blue lines
% represent the posterior expectations while the shaded areas correspond to
% 90% credible intervals.

i       = find(ismember({data.country},country)); % country index
U       = [1,4,5];
spm_COVID_ci(DCM{i}.Ep,DCM{i}.Cp,[],U);

% add seasonal flu rates
%--------------------------------------------------------------------------
subplot(2,2,2), hold on
x       = get(gca,'XLim');
plot(x,[FLU(1) FLU(1)],'-.r',x,[FLU(2) FLU(2)],'-.r')


% demonstrate routines: mitigation strategies
%==========================================================================
spm_figure('GetWin','Mitigation: posterior predictions'); clf;
%--------------------------------------------------------------------------

% get country and priors
%--------------------------------------------------------------------------
c     = find(ismember({data.country},country)); % country index

% different policies under different kinds of immunity
%--------------------------------------------------------------------------
Tim   = [4,32];                             % period of immunity (months)
tes   = [1/8,1];                            % selective testing
sde   = [1/32 1/4];                         % social distance threshold
for i = 1:numel(Tim)
    for j = 1:numel(tes)
        for k = 1:numel(sde)
            
            % set parameters
            %--------------------------------------------------------------
            P      = DCM{c}.Ep;
            P.Tim  = log(Tim(i));
            P.tes  = log(tes(j));
            P.sde  = log(sde(k));
            
            % evaluate credible interval for cumulative deaths
            %--------------------------------------------------------------
            pol{i,j,k} = sprintf('IM(%d)-ST(%d)-SD(%d)',i,j,k);
            [S,CS]     = spm_COVID_ci(P,DCM{c}.Cp,DCM{c}.Y,1);
            
            % save predictive posterior over final values
            %--------------------------------------------------------------
            SE(i,j,k)  = full(S(end));
            SC(i,j,k)  = full(CS(end,end));
            
        end
    end
end


% plot results
%--------------------------------------------------------------------------
subplot(2,1,2)
spm_plot_ci(SE(:),SC(:))
ylabel('Mitigation: cumulative deaths','FontSize',16)
set(gca,'XTick',1:numel(pol),'XTickLabel',pol)
camorbit(90,0), axis square, box off

% demonstrate routines: predictive validity
%==========================================================================
% (predictive validity - early). This figure uses the same format as above;
% however, here, the posterior estimates are based upon partial data, from
% early in the timeseries for an exemplar country. These estimates are
% based upon the empirical priors following parametric empirical Bayes. The
% red dots show the outcomes that were observed but not used to estimate
% the expected trajectories (or confidence intervals). This example
% illustrates the predictive validity of the estimates for a 10 day period
% following the last datapoint. This captures the rise to the peak of new
% cases in Italy.

% 10 day ahead forecast for Italy 
%--------------------------------------------------------------------------
i  = find(ismember({data.country},'Italy'));
spm_COVID_PV(DCM,i,10);

return


% auxiliary routines
%__________________________________________________________________________

% table of posterior estimates
%==========================================================================
[pE,pC,str] = spm_COVID_priors;
data        = DATA_COVID_JHU;

i   = find(ismember({data.country},'United Kingdom'));
Y   = DCM{i}.Y;
Ep  = DCM{i}.Ep;
Cp  = DCM{i}.Cp;
ep  = spm_vec(Ep); Cp = diag(Cp);
Tab = {};
c   = spm_invNcdf(0.05);
for i = 1:numel(str.names)
    Tab{i,1} = str.names{i};
    Tab{i,2} = exp(ep(i));
    Tab{i,3} = exp(ep(i) + c*sqrt(Cp(i)));
    Tab{i,4} = exp(ep(i) - c*sqrt(Cp(i)));
end
Table  = cell2table(Tab)
writetable(Table,'Table','FileType','spreadsheet');


% peaks
%--------------------------------------------------------------------------
i     = find(ismember({data.country},'United Kingdom'));
Z     = spm_COVID_gen(Ep,[],1:3);
for i = 1:size(Z,2)
    [d,j]  = max(Z(:,i));
    disp(j - length(Y)), disp('days to peak ')
    j      = find(Z(:,i) > 16,1,'last');
    disp(j - length(Y)), disp('days to < 16')
end


% demonstrate routines: face validation of inversion scheme
%==========================================================================

% Gaussian priors over model parameters
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% Sample model parameters from prior distribution
%--------------------------------------------------------------------------
dP    = diag(sqrt(spm_vec(pC)))*randn(spm_length(pE),1);
P     = spm_unvec(spm_vec(pE) + dP,pE);
[Y,X] = spm_COVID_gen(P,M,1:2);

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
[Z,X] = spm_COVID_gen(Ep,M,1:4);
spm_COVID_plot(Z,X,Y)


