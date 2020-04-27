function [P,C,str,rfx] = spm_COVID_priors
% Generate prior expectation and covariance log parameters
% FORMAT [pE,pC,str,rfx] = spm_COVID_priors
% 
% pE          - prior expectation (structure)
% pC          - prior covariances (structure)
% str.factor  - latent or hidden factors
% str.factors - levels of each factor
% str.outcome - outcome names (see spm_COVID_gen)
% str.names   - parameter names
% str.field   - field names of random effects
% rfx         - indices of random effects
%
% This routine assembles the (Gaussian) and priors over the parameters of a
% generative model for COVID-19. This generative model is based upon a mean
% field approximation to ensemble of population dynamics, in which four
% marginal distributions are coupled through probability transition
% matrices. The marginal distributions correspond to 4 factors;
% namely,location, infection, clinical and diagnostic or testing states.
% The parameters of this model determine the initial (probability) states
% and the transitions among the states that show certain conditional
% independences.
%
% These parameters can either be interpreted in terms of the probability of
% moving from one state to another of a given factor, conditioned on
% another. Alternatively, in some instances (specifically, staying in the
% same state),the parameters can be thought of as log transformed rate
% constants or inverse time constants.
%
% All the parameters of this generative model are log scale parameters. In
% other words, the parameters are non-negative but are  encoded in terms of
% their logarithms. This means that priors over parameters can be specified
% in terms of a prior expectation and covariance and Gaussian assumptions
% (i.e., lognormal priors over scale parameters).
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_COVID_priors.m 7841 2020-04-27 18:18:27Z karl $

% sources and background
%--------------------------------------------------------------------------
% https://jamanetwork.com/journals/jama/fullarticle/2761044
% https://www.statista.com/chart/21105/number-of-critical-care-beds-per-100000-inhabitants/
%--------------------------------------------------------------------------
% The NHS also maintains critical care beds for patients who are seriously
% ill and require constant support. Unlike most other categories of
% hospital bed, the total number of critical care beds has increased in
% recent years. In 2011/12 there were around 5,400 critical care beds, by
% 2019/20 this had risen to 5,900 (NHS England 2019b) (Figure 5). Of these,
% around 70 per cent are for use by adults and the remainder for children
% and infants.
%--------------------------------------------------------------------------
% https://www.gov.uk/guidance/coronavirus-covid-19-information-for-the-public
%
% As of 9am on 28 March 2020, a total of 120,776 people have been tested,
% of which 103,687 were confirmed negative and 17,089 were confirmed
% positive.
%--------------------------------------------------------------------------
% As of 5pm on 27 March 2020, 1,019 patients in the UK who tested positive
% for coronavirus (COVID-19) have died.
%--------------------------------------------------------------------------

% parameter names (where %** denotes fixed effects)
%==========================================================================
names{1}  = 'initial cases'; %**
names{2}  = 'size of population';
names{3}  = 'initial immunity';
names{4}  = 'P(work | home)';
names{5}  = 'social distancing';
names{6}  = 'bed availability';  
names{7}  = 'contacts: home';
names{8}  = 'contacts: work';
names{9}  = 'P(contagion | contact)';
names{10} = 'infected period';
names{11} = 'contagious period';
names{12} = 'incubation period';
names{13} = 'P(ARDS | symptoms)';
names{14} = 'symptomatic period';
names{15} = 'acute RDS period';
names{16} = 'P(fatality | CCU)';
names{17} = 'P(survival | home)';
names{18} = 'threshold for testing'; %**
names{19} = 'test rate'; %**
names{20} = 'test delay'; %**
names{21} = 'P(tested | uninfected)';  %**
names{22} = 'immune period';  %**
names{23} = 'back to work policy';  %**
names{24} = 'resistant';  %**



% random effects (i.e., effects that are common in countries)
%--------------------------------------------------------------------------
rfx       = [2:17];

% latent or hidden factors
%--------------------------------------------------------------------------
factors   = {'Location','Infection','Clinical','Testing'};

factor{1} = {'home','work','CCU','morgue'};
factor{2} = {'susceptible','infected','infectious','immune','resistant'};
factor{3} = {'none','symptoms','ARDS','deceased'};
factor{4} = {'untested','waiting','positive','negative'};

% labels or strings for plotting
%--------------------------------------------------------------------------
str.outcome = {'Death rate','New cases','CCU occupancy','ERR','Herd immunity'};
str.factors = factors;
str.factor  = factor;
str.names   = names;

% cut and paste to see the effects of changing different prior expectations
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if false
    pE    = spm_COVID_priors; M.T = 180;
    [Y,X] = spm_COVID_gen(pE,M,3); u = exp(pE.cap + pE.N)*1e6;
    spm_COVID_plot(Y,X,[],u)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% Expectations (either heuristic or taken from the above sources)
%==========================================================================
P.n   = 1;                    % number of initial cases
P.N   = 4;                    % population size (in millions)
P.m   = 1e-6;                 % herd immunity (proportion)

% location parameters
%--------------------------------------------------------------------------
P.out = 1/3;                  % P(work | home)
P.sde = 1/32;                 % social distancing threshold
P.cap = 16/100000;            % bed availability threshold (per capita)

% infection parameters
%--------------------------------------------------------------------------
P.Rin = 4;                    % effective number of contacts: home
P.Rou = 48;                   % effective number of contacts: work
P.trn = 1/4;                  % P(contagion | contact)
P.Tin = 4;                    % infected (pre-contagious) period
P.Tcn = 4;                    % infectious (contagious) period

% clinical parameters
%--------------------------------------------------------------------------
P.Tic = 8;                    % incubation period (days)
P.sev = 1/128;                % P(severe symptoms | symptomatic)
P.Tsy = 5;                    % symptomatic period
P.Trd = 12;                   % acute RDS   period
P.fat = 1/2;                  % P(fatality | CCU)
P.sur = 1/16;                 % P(survival | home)

% testing parameters
%--------------------------------------------------------------------------
P.tft = 1/1024;               % threshold: testing capacity
P.sen = 1/1024;               % rate:      testing capacity
P.del = 2;                    % delay:     testing capacity
P.tes = 1/8;                  % P(tested | uninfected)

% immunity
%--------------------------------------------------------------------------
P.Tim = 32;                   % period of immunity (months)
P.btw = 1/1024;               % back to work policy
P.r   = 1/1024;               % proportion resistant cases

% total mortality rate (for susceptible population)
%--------------------------------------------------------------------------
% (1 - (1 - P.dev)^(P.Tin + P.Tcn))*P.sev*P.fat*100


% Variances (mildly informative priors, apart from initial cases and size)
%==========================================================================
C.n   = 1;                    % number of initial cases
C.N   = 1;                    % size of population with mixing
C.m   = 0;                    % herd immunity (proportion)

% location parameters
%--------------------------------------------------------------------------
V     = 128;
C.out = 1/V;                  % P(going home | work)
C.sde = 1/V;                  % social distancing threshold
C.cap = 1/V;                  % bed availability threshold (per capita)

% infection parameters
%--------------------------------------------------------------------------
C.Rin = 1/V;                  % effective number of contacts: home
C.Rou = 1/V;                  % effective number of contacts: work
C.trn = 1/V;                  % P(transmission | infectious)
C.Tin = 1/V;                  % infected (pre-contagious) period
C.Tcn = 1/V;                  % contagious period

% clinical parameters
%--------------------------------------------------------------------------
C.dev = 1/V;                  % P(developing symptoms | infected)
C.sev = 1/V;                  % P(severe symptoms | symptomatic)
C.Tsy = 1/V;                  % symptomatic period
C.Trd = 1/V;                  % acute RDS   period
C.fat = 1/V;                  % P(fatality | CCU)
C.sur = 1/V;                  % P(fatality | home)

% testing parameters
%--------------------------------------------------------------------------
C.tft = 1/4;                  % threshold for testing capacity
C.sen = 1/4;                  % sensitivity:  testing capacity
C.del = 1/4;                  % delay:        testing capacity
C.tes = 1/4;                  % P(testing | uninfected)

% immunity
%--------------------------------------------------------------------------
C.Tim = 1/V;                  % period of immunity
C.btw = 0;                    % back to work policy
C.r   = 0;                    % proportion of people not susceptible



% log transform
%==========================================================================
P         = spm_vecfun(P,@log);

% field names of random effects
%--------------------------------------------------------------------------
field     = fieldnames(P);
str.field = field(rfx);

return