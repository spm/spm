function [P,C,str] = spm_SARS_priors
% Generate prior expectation and covariance log parameters
% FORMAT [pE,pC,str,rfx] = spm_SARS_priors
% 
% pE          - prior expectation (structure)
% pC          - prior covariances (structure)
% str.factor  - latent or hidden factors
% str.factors - levels of each factor
% str.outcome - outcome names (see spm_SARS_gen)
% str.names   - parameter names
% str.field   - field names of random effects
% rfx         - indices of random effects
%
% This routine assembles the (Gaussian) and priors over the parameters of a
% generative model for SARS-19. This generative model is based upon a mean
% field approximation to ensemble of population dynamics, in which four
% marginal distributions are coupled through probability transition
% matrices. The marginal distributions correspond to 4 factors; namely,
% location, infection, symptom and testing (LIST) states. The parameters of
% this model determine the initial (probability) states and the transitions
% among the states that show certain conditional independences.
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
% $Id: spm_SARS_priors.m 7939 2020-09-09 11:02:14Z karl $

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
% for coronavirus (SARS-19) have died.
%--------------------------------------------------------------------------

% parameter names (where %** denotes fixed effects)
%==========================================================================
names{1}  = 'population size';
names{2}  = 'initial cases';
names{3}  = 'unexposed';
names{4}  = 'pre-existing immunity';
names{5}  = 'initial exposed';

% location (exposure) parameters
%--------------------------------------------------------------------------
names{6}  = 'P(leaving home)';
names{7}  = 'threshold: lockdown';  
names{8}  = 'threshold: containment';
names{9}  = 'viral spreading';
names{10} = 'bed availability';
names{11} = 'sensitivity';

% infection (transmission) parameters
%--------------------------------------------------------------------------
names{12} = 'contacts: home';
names{13} = 'contacts: work';
names{14} = 'contacts: rural';
names{15} = 'transmission strength';
names{16} = 'infected period';
names{17} = 'infectious period';
names{18} = 'seropositive immunity';
names{19} = 'seronegative immunity';

% clinical parameters
%--------------------------------------------------------------------------
names{20} = 'incubation period';
names{21} = 'symptomatic period';
names{22} = 'hospitalisation period';
names{23} = 'P(ARDS|symptoms)';
names{24} = 'P(fatality|CCU): early';
names{25} = 'P(fatality|CCU): late';
names{26} = 'P(survival|home)';

% testing parameters
%--------------------------------------------------------------------------
names{27} = 'FTTI efficacy';
names{28} = 'testing: seroprevalence';
names{29} = 'testing: capacity';
names{30} = 'test selectivity';
names{31} = 'test delay';

names{32} = 'buildup testing';
names{33} = 'buildup latency';
names{34} = 'buildup period';


% latent or hidden factors
%--------------------------------------------------------------------------
factors   = {'Location','Infection','Symptoms','Testing'};

factor{1} = {'home','work','hospital','removed','isolated'};
factor{2} = {'susceptible','infected','infectious','AB +ve','AB -ve'};
factor{3} = {'none','symptoms','severe','deceased'};
factor{4} = {'untested','waiting','PCR +ve','PCR -ve'};

% labels or strings for plotting
%--------------------------------------------------------------------------
str.outcome = {'Death rate',...
               'New tests',...
               'CCU occupancy',...
               'ERR',...
               'Herd immunity',...
               'Tests',...
               'Contagion risk (%)',...
               'Prevalence {%}',...
               'FTTI target',...
               'New cases'};
str.factors = factors;
str.factor  = factor;
str.names   = names;

% cut and paste to see the effects of changing different prior expectations
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if false
    pE    = spm_SARS_priors; M.T = 365;
    [Y,X] = spm_SARS_gen(pE,M,1:3);
    u     = exp(pE.cap + pE.N)*1e6;
    spm_SARS_plot(Y,X,[],u)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% Expectations (either heuristic or taken from the above sources)
%==========================================================================
P.N   = 64;                   % (01) population size (millions)
P.n   = 64;                   % (02) initial cases (cases)
P.m   = 1/2;                  % (03) unexposed (proportion)
P.r   = 1/2;                  % (04) pre-existing immunity (proportion)
P.o   = 1/16;                 % (05) initial exposed proportion

% location (exposure) parameters
%--------------------------------------------------------------------------
P.out = 1/3;                  % (06) P(leaving home)
P.sde = 1/64;                 % (07) lockdown threshold
P.qua = 1/256;                % (08) Quarantine threshold
P.exp = 32;                   % (09) contacts (viral spreading)
P.cap = 8/100000;             % (10) bed availability (per capita)
P.s   = 4;                    % (11) sensitivity parameter

% infection (transmission) parameters
%--------------------------------------------------------------------------
P.Nin = 5;                    % (12) effective number of contacts: home
P.Nou = 48;                   % (13) effective number of contacts: work
P.Nru = 32;                   % (14) effective number of contacts: rural
P.trn = 1/3;                  % (15) transmission strength
P.Tin = 4;                    % (16) infected period (days)
P.Tcn = 4;                    % (17) infectious period (days)
P.Tim = 32;                   % (18) seropositive immunity (months)
P.res = .4;                   % (19) seronegative immunity (proportion)

% clinical parameters
%--------------------------------------------------------------------------
P.Tic = 12;                   % (20) incubation period (days)
P.Tsy = 8;                    % (21) symptomatic period (days)
P.Trd = 12;                   % (22) CCU period (days)
P.sev = 1/32;                 % (23) P(ARDS | symptomatic)
P.fat = 1/2;                  % (24) P(fatality | ARDS, CCU): early
P.lat = 1/4;                  % (25) P(fatality | ARDS, CCU): late
P.sur = 1/8;                  % (26) P(survival | ARDS, home)

% testing parameters
%--------------------------------------------------------------------------
P.ttt = 1/16;                 % (27) FTTI efficacy
P.tts = 1;                    % (28) testing: seroprevalence
P.lim = 2/1000;               % (29) testing: capacity
P.tes = 4;                    % (30) test selectivity (for infection)
P.del = 4;                    % (31) test delay (days)

P.sus = 2/1000;               % (32) buildup testing
P.ont = 2;                    % (33) buildup latency (months)
P.stt = 8;                    % (34) buildup period (days)

% infection fatality (for susceptible population)
%--------------------------------------------------------------------------
% IFR (hospital): P.sev*P.fat*100
% IFR (carehome): P.sev*(1 - P.sur)*100
% IFR (hospital): exp(Ep.sev)*exp(Ep.fat)*100
% IFR (carehome): exp(Ep.sev)*(1 - exp(Ep.sur))*100

% Variances (mildly informative priors, apart from initial cases and size)
%==========================================================================
U     = 1;                    % flat priors
V     = 1/16;                 % informative priors
W     = 1/256;                % precise priors

C.N   = U;                    % (01) population size (millions)
C.n   = U;                    % (02) initial cases (cases)
C.m   = W;                    % (03) unexposed (proportion)
C.r   = W;                    % (04) pre-existing immunity (proportion)
C.o   = W;                    % (05) initial exposed proportion


% location (exposure) parameters
%--------------------------------------------------------------------------
C.out = W;                    % (06) P(leaving home)
C.sde = W;                    % (07) lockdown threshold
C.qua = W;                    % (08) Quarantine threshold
C.exp = W;                    % (09) contacts (viral spreading)
C.cap = W;                    % (10) bed availability (per capita)
C.s   = W;                    % (11) sensitivity parameter

% infection (transmission) parameters
%--------------------------------------------------------------------------
C.Nin = V;                    % (12) effective number of contacts: home
C.Nou = W;                    % (13) effective number of contacts: work
C.Nru = W;                    % (14) effective number of contacts: rural
C.trn = W;                    % (15) transmission strength
C.Tin = W;                    % (16) infected period (days)
C.Tcn = W;                    % (17) infectious period (days)
C.Tim = W;                    % (18) seropositive immunity (months)
C.res = W;                    % (19) seronegative immunity (proportion)

% clinical parameters
%--------------------------------------------------------------------------
C.Tic = W;                    % (20) incubation period (days)
C.Tsy = W;                    % (21) symptomatic period (days)
C.Trd = W;                    % (22) CCU period (days)
C.sev = W;                    % (23) P(ARDS | symptomatic)
C.fat = W;                    % (24) P(fatality | ARDS, CCU): early
C.lat = W;                    % (25) P(fatality | ARDS, CCU): late
C.sur = W;                    % (26) P(survival | ARDS, home)

% testing parameters
%--------------------------------------------------------------------------
C.ttt = V;                    % (27) FTTI efficacy
C.tts = W;                    % (28) testing: seroprevalence
C.lim = W;                    % (29) testing capacity
C.tes = W;                    % (30) test selectivity (for infection)
C.del = W;                    % (31) test delay (days)

C.sus = W;                    % (32) buildup testing
C.ont = U;                    % (33) buildup latency (months)
C.stt = U;                    % (34) buildup period (days)


% log transform
%==========================================================================
P         = spm_vecfun(P,@log);

% field names of random effects
%--------------------------------------------------------------------------
str.field = fieldnames(P);

return