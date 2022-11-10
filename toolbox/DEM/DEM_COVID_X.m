function [DCM] = DEM_COVID_X(Y,Data)
% FORMAT [DCM] = DEM_COVID_X(Y,data)
% data    - data to model [default: data = DATA_COVID_US]
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% This routine illustrates the Bayesian model inversion of a generative
% model of coronavirus spread using variational techniques (variational
% Laplace). This (pandemic) model is composed of regional (epidemic)
% models. In brief, the model for a single region comprises four factors,
% each with four states, giving 256 states or compartments per region.
% These regional models are then assembled to model the coupling among
% eight regions giving 256^8 compartments. However, due to certain
% conditional independencies, this can be treated as a collection of 256
% compartmental models; providing one carefully links the state of one
% region to the state of another. Here, this linking or connectivity is
% parameterised in terms of a probability flux or exchange of people from
% one regional population to another. Regional factors include location,
% immune status, clinical status and testing status. The transitions among
% states of any factor depends upon other factors. For example, the
% probability that I will move from a state of being asymptomatic to being
% symptomatic depends upon whether I am infected or not. Similarly, the
% probability that I will move from one region to another depends upon
% whether I am at work (i.e., not at home). In short, the exchange between
% different regional populations is limited to the people who are not at
% home and are consequently in a position to travel. The parameters of
% interregional coupling correspond to rate constants or effective
% connectivity that can be reciprocal and asymmetric. For example, the
% probability of moving to New York from New Jersey does not have to be the
% same as a probability of moving from New Jersey to New York. Note that
% the movement between regions can be restricted to a chain. In other
% words, to get from the first state to the last state, I have to go
% through all other states.
%
% Each subsection produces one or two figures that are described in the
% annotated (Matlab) code. These subsections call various subroutines that
% provide a more detailed description of things like the generative model,
% its priors and the evaluation of confidence intervals.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% Get data (see DATA_COVID_US): an array with a structure for each State
%==========================================================================
if nargin < 1, [Y,Data] = DATA_COVID_US; end

% number of regions (e.g.,US states)
%--------------------------------------------------------------------------
i    = 1:10;                             % number of regions
Y    = Y(:,:,i);                         % empirical cases
data = Data(i);                          % Meta data

% Inversion (i.e., fitting) of empirical data for each region
%==========================================================================
Fsi   = spm_figure('GetWin','SI'); clf;

% assemble (Gaussian) priors over model parameters
%----------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% Bayesian inversion (placing posteriors in a cell array of structures)
%----------------------------------------------------------------------
for i = 1:numel(data)
    
    % variational Laplace (estimating log evidence (F) and posteriors)
    %------------------------------------------------------------------
    set(Fsi,'name',data(i).state)
    pE.N      = log(data(i).pop*1e-6);   % fix (log) population size
    pC.N      = 0;                       % with no uncertainty  
    [F,Ep,Cp] = spm_COVID(Y(:,:,i),pE,pC);
    
    % assemble prior and posterior estimates (and log evidence)
    %------------------------------------------------------------------
    M.Q(i).Ep = Ep;
    M.Q(i).Cp = Cp;
    M.Q(i).F  = F;
    
end

% assemble interregional but priors
%======================================================================
[pE,pC,STR,erc] = spm_COVID_priors_R(data);
str.names       = STR.names;
str.regions     = STR.regions;

% complete model specification
%----------------------------------------------------------------------
M.G    = @spm_COVID_US;            % generative function
M.FS   = @(Y)sqrt(Y);              % feature selection  (link function)
M.pE   = pE;                       % prior expectations (parameters)
M.pC   = pC;                       % prior covariances  (parameters)
M.hE   = 0;                        % prior expectation  (log-precision)
M.hC   = 1/64;                     % prior covariances  (log-precision)
M.T    = size(Y,1);                % number of samples  (time)
M.data = data;                     % number of samples  (regions)
M.erc  = erc;                      % effective regional connectivity
U      = 1:size(Y,2);              % number of response variables

% model inversion with Variational Laplace (Gauss Newton) {1}
%======================================================================
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);

% assemble prior and posterior estimates (and log evidence)
%----------------------------------------------------------------------
DCM.M  = M;
DCM.Ep = Ep;
DCM.Cp = Cp;
DCM.Eh = Eh;
DCM.F  = F;
DCM.Y  = Y;

% save
%--------------------------------------------------------------------------
clear ans Fsi, save COVID_X

% region specific predictions
%==========================================================================                                % get estimated parameters
% This figure reports the predicted outcomes and underlying latent states
% generating those outcomes over a one-year period. The upper panels show
% some key outcomes, some of which are measurable. Here, daily death rates
% are shown in blue, new cases in red and CCU occupancy in orange. The
% lines correspond to the predictions of the model, while the dots are
% empirical data available at the time of fitting. The dotted line in the
% upper right panel corresponds to the typical critical care capacity of a
% large city. The same results are shown on the upper right panel in terms
% of cumulative new cases (red) and deaths (blue). The underlying or latent
% causes of this mortality are shown in the lower panels. These are
% organised according to the four factors of the generative (dynamic
% causal) model. In each panel, the latent states are plotted for the
% regions considered in this analysis. The location factor shows that,
% under this strategy, the number of people away from the home (or
% equivalent location) decreases sharply at the onset of the outbreak and
% then recovers slowly over the ensuing weeks. In terms of infection, there
% is a rapid acquisition of immunity over the first months of the pandemic.
% At any one time, about 10% or less of the population is either infected
% or infectious. In terms of the clinical expression of these infections,
% 10% or less of people will experience symptoms and a small minority will
% progress to acute respiratory distress, from which they may recover or
% die. Under this model, positive test results for the virus (based on
% buccal swabs) accumulate over time as more and more people are tested. In
% the initial phases of the outbreak, most people are negative. However,
% during the onset of the pandemic about a half to a third of people tested
% are positive. This proportion declines over the ensuing months.
%--------------------------------------------------------------------------
DCM.M.T = 365;
spm_figure('GetWin','USA - 32 months immunity'); clf;

[Y,X]   = spm_COVID_US(DCM.Ep,DCM.M,1:3);
spm_COVID_plot(Y,X,DCM.Y)

% repeat with the loss of immunity
%--------------------------------------------------------------------------
spm_figure('GetWin','USA - 3 months immunity'); clf;
%--------------------------------------------------------------------------
T     = DCM.M;
for k = 1:numel(T.Q)
    T.Q(k).Ep.Tim = log(4);
end

[Y,X] = spm_COVID_US(DCM.Ep,T,1:3);
spm_COVID_plot(Y,X,DCM.Y)

% Differences among countries, in terms of parameters
%==========================================================================
spm_figure('GetWin','Within region parameters'); clf;
%--------------------------------------------------------------------------
% (differences among regions). This figure reports the differences among
% regions in terms of selected parameters of the generative model,
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
%--------------------------------------------------------------------------
[PE,PC,Pnm] = spm_COVID_priors;
% names{1}  = 'initial cases'; %**
% names{2}  = 'size of population';
% names{3}  = 'initial immunity';
% names{4}  = 'P(work | home)';
% names{5}  = 'social distancing';
% names{6}  = 'bed availability';  
% names{7}  = 'contacts: home';
% names{8}  = 'contacts: work';
% names{9}  = 'P(contagion | contact)';
% names{10} = 'infected period';
% names{11} = 'contagious period';
% names{12} = 'incubation period';
% names{13} = 'P(ARDS | symptoms)';
% names{14} = 'symptomatic period';
% names{15} = 'CCU period';
% names{16} = 'P(fatality | CCU)';
% names{17} = 'P(survival | home)';
% names{18} = 'trace and test'; %**
% names{19} = 'response testing'; %**
% names{20} = 'test delay'; %**
% names{21} = 'test selectivity'; %**
% names{22} = 'sustained testing'; %**
% names{23} = 'immune period'; %**
% names{24} = 'resistant'; %**
% names{25} = 'initial testing'; %**

% assemble parameters
%--------------------------------------------------------------------------
P     = [];                                   % posterior expectations
C     = [];                                   % posterior variances
for i = 1:numel(DCM.M.Q)
    P(:,i) = spm_vec(DCM.M.Q(i).Ep);
    C(:,i) = diag(DCM.M.Q(i).Cp);
end

% report selected parameters (see spm_COVID_priors)
%--------------------------------------------------------------------------
p     = [2,4,5,7,8,9,11,12,13,15,16,19];
for i = 1:length(p)
    
    % posterior density
    %----------------------------------------------------------------------
    subplot(4,3,i)
    Ep   = P(p(i),:);
    Cp   = C(p(i),:);
    spm_plot_ci(Ep',Cp,[],[],'exp'), hold on
    title(Pnm.names{p(i)},'FontSize',16)
    xlabel('country'), axis square, box off
    
    % region with greatest (and least) map estimate
    %----------------------------------------------------------------------
    [d,j] = max(exp(Ep));
    text(j,d,data(j).state,'FontSize',8);
    [d,j] = min(exp(Ep));
    text(j,d,data(j).state,'FontSize',8);
    
end


% Between country analysis (hierarchical or parametric empirical Bayes)
%==========================================================================
spm_figure('GetWin','BMR - all'); clf;
%--------------------------------------------------------------------------
% (connectivity and viral spread) This figure reports the connectivity
% among regions under a model of exchange between states. The upper panel
% shows the maximum a posteriori (MAP) estimates as an adjacency matrix
% that quantifies the rate or probability that a region in each column will
% deliver a proportion of its (out-of-home) population to a region in the
% rows. The corresponding log parameters for the (non-zero) connectivity
% are shown in the middle panel in terms of their posterior expectations
% (blue bars) and 90% Bayesian credible intervals (pink bars). The lower
% panel shows the posterior probability of a model with and without each of
% these parameters, based upon Bayesian model comparison.
%--------------------------------------------------------------------------
DCM.beta  = -16;
DCM.gamma = exp(-8);
[BMA,BMR] = spm_dcm_bmr_all(DCM,'erc');

subplot(3,2,2)
ERC     = exp(BMA.Ep.erc);
ERC     = ERC + ERC';
imagesc(1 - ERC)
set(gca,'YTick',1:numel(str.regions),'Yticklabel',str.regions,'FontSize',10)
xlabel('State'), title('Regional connectivity','FontSize',16),
axis image, box off

% Sensitivity analysis: which factors determine cumulative deaths?
%==========================================================================
% (sensitivity analysis) This figure reports the effect of changing each
% parameter on cumulative deaths over an 18-month period. The upper panel
% shows the rate of increase (or decrease) in cumulative deaths per unit
% change in the (log) parameters. These sensitivity metrics are based upon
% a first order Taylor expansion about the maximum a posteriori values
% shown in the lower panel. The blue bars correspond to the most likely
% estimate and pink bars report (a large number assumption approximation
% to) the 90% credible intervals.
%--------------------------------------------------------------------------
spm_COVID_R_cii(DCM,1,'USA');

% illustrate difference between hard and soft distancing rules
%==========================================================================
spm_figure('GetWin','Formal strategies'); clf;
%--------------------------------------------------------------------------
% (soft and hard strategies) This figure shows the different kinds of
% social distancing strategies that could be adopted, i.e., considered in
% the modelling. Both parameterise the degree of social distancing as a
% function of the proportion of the population that are currently infected
% (or in critical care). This proportion can, in principle, be estimated
% directly or indirectly, given current testing capabilities. Both are
% decreasing functions of the prevalence of infection. In other words, as
% the prevalence of infection increases, the probability of leaving home
% (i.e., self-isolation) decreases. The strategies are distinguished by the
% form of this decrease. The hard (threshold) strategy is based upon a
% threshold, afforded by the sigmoid function in magenta. Conversely, the
% soft (power) strategy decreases smoothly as a power function of
% prevalence.
%--------------------------------------------------------------------------
spm_sigma = @(x,u)spm_phi(4*(u - x)/u);
p         = linspace(0,32,128);
sde       = 16;

% threshold strategies
%--------------------------------------------------------------------------
subplot(3,1,1)
plot(p,spm_sigma(p/100,1/sde),p,(1 - p/100).^sde,':',[0,0] + 100/sde,[0,1],'-.');
title('Social distancing rules','FontSize',16), box off
xlabel('percentage of the population affected (%)')
ylabel('degree of social distancing (%)')
axis square, legend({'Sigmoid function','Power function'}),legend('boxoff')


% Federal versus local policy: which factors determine cumulative deaths?
%==========================================================================
spm_figure('GetWin','strategy analysis I'); clf
%--------------------------------------------------------------------------
% how long is a lockdown? This figure shows the results of simulations
% under different levels of social distancing. This illustrates the rate of
% new cases and deaths per day (in the upper left panel) and the underlying
% or latent causes (in the lower four panels). Here, the lines report the
% average rates and probabilities over regions, for different levels of
% social distancing. These levels are evaluated in 16 steps from high to
% low social distancing (i.e., low to high pressure). The subsequent figure
% summarises these results in terms of cumulative deaths, cumulative
% working days and weeks lost to lockdown, as a function of social
% distancing.
%--------------------------------------------------------------------------
M     = DCM.M;                              % model
P     = DCM.Ep;                             % expansion point
M.T   = 18*32;                              % over 18 months
sde   = linspace(-4,4,16);                  % range of social distancing
S     = sde;
for i = 1:numel(sde)
    
    % social distancing threshold
    %----------------------------------------------------------------------
    for j = 1:numel(M.Q)
        M.Q(j).Ep.sde = DCM.M.Q(j).Ep.sde + sde(i);
    end
    [Y,X] = spm_COVID_US(P,M,1:4);
    
    % average over regions
    %----------------------------------------------------------------------
    Y     = sum(Y,3);
    XX    = cell(size(X,1),1);
    for j = 1:size(X,1)
        XX{j} = 0;
        for k = 1:size(X,2)
            XX{j} = XX{j} + X{j,k}/size(X,2);
        end
    end
    Deaths(i)   = sum(Y(:,1));
    Days(i)     = sum(Y(:,4));
    Duration(i) = sum(XX{1}(:,2) < 16/100)/7;
    
    % plot results and hold graph
    %----------------------------------------------------------------------
    if rem(i,2)
        spm_COVID_plot(Y(:,1:3),XX)
        for j = 1:6
            subplot(3,2,j), hold on
            set(gca,'ColorOrderIndex',1);
        end
    end
end

% cumulative deaths, working days and lost weeks
%--------------------------------------------------------------------------
spm_figure('GetWin','strategy analysis II'); clf

subplot(2,2,1), hold off
plot(sde,Deaths,':',sde,Deaths,'o',[0,0],[min(Deaths),max(Deaths)],'-.')
title('Cumulative deaths','FontSize',16),
xlabel('social distancing threshold')
ylabel('cumulative deaths')
axis square, box off

subplot(2,2,2), hold off
plot(sde,Days,':',sde,Days,'o',[0,0],[min(Days),max(Days)],'-.')
title('Working days','FontSize',16),
xlabel('social distancing threshold')
ylabel('cumulative deaths')
axis square, box off

subplot(2,2,3), hold off
plot(sde,Duration,':',sde,Duration,'o',[0,0],[min(Duration),max(Duration)],'-.')
title('Lost weeks','FontSize',16),
xlabel('social distancing threshold')
ylabel('weeks')
axis square, box off

subplot(2,2,4), hold off
i = Duration(find(sde > 0,1,'first'));
plot(Duration,Deaths,':',Duration,Deaths,'o',[i,i],[min(Deaths),max(Deaths)],'-.')
title('Exit strategies','FontSize',16),
xlabel('weeks lost')
ylabel('cumulative deaths')
axis square, box off

fprintf('%0.2f weeks\n',i)


%==========================================================================
spm_figure('GetWin','Policies'); clf; clear P S
%--------------------------------------------------------------------------
% (different strategies evaluated) These bar charts report the effects of
% different strategies on cumulative deaths (left column), total number of
% working days (middle column) and peak occupancy of CCU (right column).
% The top row shows the results based upon the MAP parameters. The middle
% row reproduces the same analysis but under more pessimistic assumptions
% about the retention of immunity. Specifically, we reduce the time
% constant for retaining immunity from 32 months to 4 months. The lower row
% shows the equivalent results when including a back to work policy based
% upon serological testing. The only effect of this, under the current
% model, is to reverse the effect of a hard versus partial social
% distancing strategy on the number of working days lost. This follows
% because returning people who are immune to the community has no effect on
% morbidity ordinand for critical care; however, it does take the pressure
% off the economy.
%--------------------------------------------------------------------------
DCM.M.T = 18*32;                        % over 18 months
out     = {'Deaths','Days','Peak CCU'};
policy  = {'lockdown','relaxed'};


for q = 1:3
    
    % different kinds of immunity
    %----------------------------------------------------------------------
    M     = DCM.M;                      % model (within region parameters)
    Ep    = DCM.Ep;                     % between region parameters
    sde   = [ 0  2];                    % social distance threshold
    fed   = [-8 -1/8];                  % federal policy
    if q == 2
        for k = 1:numel(M.Q)
            M.Q(k).Ep.Tim = log(4);
        end
    end
    if q == 3
        for k = 1:numel(M.Q)
            M.Q(k).Ep.btw = -1/8;
        end
    end
    
    % different policies under different kinds of immunity
    %----------------------------------------------------------------------
    for i = 1:numel(sde)
        for j = 1:numel(fed)
            
            % set parameters
            %==============================================================
            
            % social distancing threshold
            %--------------------------------------------------------------
            for k = 1:numel(M.Q)
                M.Q(k).Ep.sde = DCM.M.Q(k).Ep.sde + sde(i);
            end
            
            % Federal parameters
            %--------------------------------------------------------------
            Ep.fed = fed(j);
            
            % simulate pandemic and average over regions
            %--------------------------------------------------------------
            Y      = spm_COVID_US(Ep,M,1:4);
            Y      = sum(Y,3);
            
            % outcome measures
            %--------------------------------------------------------------
            S(i,j,1) = sum(Y(:,1));            % cumulative deaths
            S(i,j,2) = sum(Y(:,4));            % cumulative working days
            S(i,j,3) = max(Y(:,3));            % peak CCU occupancy
            
        end
    end
    
    % cumulative deaths cases, working days and CCU occupancy
    %----------------------------------------------------------------------
    for i = 1:size(S,3)
        subplot(3,3,i + (q - 1)*3), bar(S(:,:,i)), hold on
        plot([0,2] + 1/2,[0,0] + max(max(S(:,:,i))),'-.')
        set(gca,'XTick',1:numel(policy),'Xticklabel',policy)
        ylabel('incidence'), title(out{i},'FontSize',16)
        xlabel('social distancing threshold')
        axis square, box off
    end
    if q > 2, legend({'State policy','Federal policy'}), legend('boxoff'), end
    
end



% prevalence of immunity in California
%==========================================================================
i  = ismember([str.regions],'California');
if ~any(i), return, end

spm_figure('GetWin','CA'); clf;
%--------------------------------------------------------------------------
% We measured the seroprevalence of antibodies to SARS-CoV-2 in Santa Clara
% County. Methods On 4/3-4/4, 2020, we tested county residents for
% antibodies to SARS-CoV-2 using a lateral flow immunoassay.The unadjusted
% prevalence of antibodies to SARS-CoV-2 in Santa Clara County was 1.5%
% (exact binomial 95CI 1.11-1.97%), and the population-weighted prevalence
% was 2.81% (95CI 2.24-3.37%). Under the three scenarios for test
% performance characteristics, the population prevalence of COVID-19 in
% Santa Clara ranged from 2.49% (95CI 1.80-3.17%) to 4.16% (2.58-5.70%)
%--------------------------------------------------------------------------
T     = datenum('03-Apr-2020') - datenum('20-Jan-2020') - 16;
i     = ismember([str.regions],'California');
[Y,X] = spm_COVID_US(DCM.Ep,DCM.M,5);
X     = X{2,i};

% get infection status predictions for this state
%--------------------------------------------------------------------------
j     = size(DCM.Y,1);
X     = X(1:j,2:4);

subplot(2,1,1)
plot(X*100)
ylabel('prevalence (%)'), title('infection status','FontSize',16)
xlabel('time (days)'), set(gca,'YLim',[0 16])
axis square, box off

% superimpose empirical estimates
%--------------------------------------------------------------------------
t   = T; hold on
plot([t + 1,t + 1],[1.80,3.17],[t,t],[2.58,5.70],'Linewidth',8)
legend({'infected','infectious','immune','prevalence','prevalence'})


%--------------------------------------------------------------------------
% The initial results from the first large-scale study tracking the spread
% of the coronavirus in the county found that 4.1% of adults have
% antibodies to the virus in their blood, an indication of past exposure.
% That translates to roughly 221,000 to 442,000 adults who have recovered
% from an infection, once margin of error is taken into account, according
% to the researchers conducting the study. The county had reported fewer
% than 8,000 cases at that time.
% The findings suggest the fatality rate may be much lower than previously
% thought. But although the virus may be more widespread, the infection
% rate still falls far short of herd immunity that, absent a vaccine, would
% be key to return to normal life
% https://www.latimes.com/california/story/2020-04-20/coronavirus-serology-testing-la-county

cases = cumsum(DCM.Y(:,2,i));
j     = find(cases < 8000,1,'last');
plot(j,4.1,'.r','MarkerSize',32)
legend({'infected','infectious','immune','prevalence','prevalence','LA Times'})


return



% auxiliary functions
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


% Second waves in China
%==========================================================================
spm_figure('GetWin','China'); clf;
%--------------------------------------------------------------------------
% (second waves in China) This figure illustrates a second wave of new cases
% in reports from China. The docs represent empirical records of new cases
% and deaths as a function of weeks from the onset of the outbreak. The
% lines correspond to the predicted incidences, under a (single region)
% model described in DEM_COVID. These data are presented to illustrate
% a resurgence of new cases several weeks after the first wave that,
% crucially, is not accompanied by a secondary increase in death rates.
%--------------------------------------------------------------------------
i       = find(ismember({Data.country},'China'));          % country index
M.T     = 100;
[Y,X]   = spm_COVID_gen(DCM{i}.Ep,M,1:2);
spm_COVID_plot(Y,X,DCM{i}.Y);
subplot(3,2,1), set(gca,'YLim',[0 800])

return

% Predictive validity for London: T-day lookahead projection
%==========================================================================
spm_figure('GetWin','UK'); clf;
%--------------------------------------------------------------------------
load COVID_DCM
Data    = DATA_COVID_JHU;

% predictive validity for London (UK)
%--------------------------------------------------------------------------
i       = find(ismember({Data.country},'United Kingdom')); % country index
M.T     = 180;                                             % six-month period
T       = 10;                                              % 4-Apr-20
DCM{i}.M.pE.N = 0;
spm_COVID_PV(DCM,i,T);

return
