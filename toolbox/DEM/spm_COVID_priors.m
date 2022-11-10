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
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

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
names{2}  = 'population size';
names{3}  = 'initial proportion';
names{4}  = 'P(work | home)';
names{5}  = 'social distancing';
names{6}  = 'bed availability';  
names{7}  = 'contacts: home';
names{8}  = 'contacts: work';
names{9}  = 'transmission strength';
names{10} = 'infected period';
names{11} = 'contagious period';
names{12} = 'incubation period';
names{13} = 'P(ARDS | symptoms)';
names{14} = 'symptomatic period';
names{15} = 'time in CCU';
names{16} = 'P(fatality | CCU)';
names{17} = 'P(survival | home)';
names{18} = 'track and trace'; %**
names{19} = 'testing latency'; %**
names{20} = 'test delay'; %**
names{21} = 'test selectivity'; %**
names{22} = 'sustained testing'; %**
names{23} = 'baseline testing'; %**
names{24} = 'immune period'; %**
names{25} = 'exempt period'; %**
names{26} = 'resistance'; %**
names{27} = 'innate immunity'; %**

% random effects (i.e., effects that are common in countries)
%--------------------------------------------------------------------------
rfx       = 2:17;

% latent or hidden factors
%--------------------------------------------------------------------------
factors   = {'Location','Infection','Symptoms','Testing'};

factor{1} = {'home','work','hospital','removed','isolated'};
factor{2} = {'susceptible','infected','infectious','Ab +ve','Ab -ve'};
factor{3} = {'none','symptoms','ARDS','deceased'};
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
    pE    = spm_COVID_priors; M.T = 365;
    [Y,X] = spm_COVID_gen(pE,M,1:3);
    u     = exp(pE.cap + pE.N)*1e6;
    spm_COVID_plot(Y,X,[],u)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% https://annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported

% Expectations (either heuristic or taken from the above sources)
%==========================================================================
P.n   = 4;                    % number of initial cases
P.N   = 8;                    % population size (in millions)
P.m   = 1/4;                  % initial proportion

% location parameters
%--------------------------------------------------------------------------
P.out = 1/3;                  % P(work | home)
P.sde = 1/32;                 % social distancing threshold
P.cap = 16/100000;            % bed availability threshold (per capita)

% infection parameters
%--------------------------------------------------------------------------
P.Rin = 4;                    % effective number of contacts: home
P.Rou = 48;                   % effective number of contacts: work
P.trn = 1/3;                  % P(contagion | contact)
P.Tin = 4;                    % infected (pre-contagious) period
P.Tcn = 4;                    % infectious (contagious) period

% clinical parameters
%--------------------------------------------------------------------------
P.Tic = 16;                   % time till symptoms (days)
P.sev = 1/32;                 % P(severe symptoms | symptomatic)
P.Tsy = 8;                    % symptomatic period
P.Trd = 10;                   % CCU period
P.fat = 1/2;                  % P(fatality | severe, CCU)
P.sur = 1/8;                  % P(survival | severe, home)

% testing parameters
%--------------------------------------------------------------------------
P.ttt = 1/10000;              % test, track and trace
P.ont = 2;                    % testing latency (months)
P.del = 4;                    % test delay (days)
P.tes = 8;                    % test selectivity (for infection)
P.sus = 4/10000;              % sustained testing
P.bas = 4/10000;              % baseline testing

% immunity
%--------------------------------------------------------------------------
P.Tim = 16;                   % period of immunity (months)
P.Tex = 2;                    % period of exemption (days)
P.r   = 1/2;                  % proportion resistant cases
P.res = 1/2;                  % proportion with innate immunity
P.tts = 16;                   % testing buildup

% total mortality rate (for susceptible population)
%--------------------------------------------------------------------------
% IFR (hospital): P.sev*P.fat*100
% IFR (carehome): P.sev*(1 - P.sur)*100
% IFR (hospital): exp(Ep.sev)*exp(Ep.fat)*100
% IFR (carehome): exp(Ep.sev)*(1 - exp(Ep.sur))*100

% We describe what we believe is the first instance of complete COVID-19
% testing of all passengers and crew on an isolated cruise ship during the
% current COVID-19 pandemic. Of the 217 passengers and crew on board, 128
% tested positive for COVID-19 on reverse transcription-PCR (59 pc). Of the
% COVID- 19-positive patients, 19 pc (24) were symptomatic; 6.2 pc (8)
% required medical evacuation; 3.1 pc (4) were intubated and ventilated;
% and the mortality was 0.8 pc (1). The majority of COVID-19-positive
% patients were asymptomatic (81 pc, 104 patients).

% probability of developing symptoms if infected (19%)
%--------------------------------------------------------------------------
% 1 - exp(-1/P.Tic)^(P.Tin + P.Tcn)

% probability of ARDS if symptomatic (4/24)
%--------------------------------------------------------------------------
% P.sev

% probability of dying if ARDS (1/4)
%--------------------------------------------------------------------------
% P.fat

% Variances (mildly informative priors, apart from initial cases and size)
%==========================================================================
U     = 1;                    % flat priors
V     = 1/16;                 % informative priors
W     = 1/256;                % precise priors

C.n   = U;                    % number of initial cases
C.N   = U;                    % size of population with mixing
C.m   = W;                    % initial proportion

% location parameters
%--------------------------------------------------------------------------
C.out = W;                    % P(going home | work)
C.sde = W;                    % social distancing threshold
C.cap = W;                    % bed availability threshold (per capita)

% infection parameters
%--------------------------------------------------------------------------
C.Rin = V;                    % effective number of contacts: home
C.Rou = V;                    % effective number of contacts: work
C.trn = V;                    % P(transmission | infectious)
C.Tin = W;                    % infected (pre-contagious) period
C.Tcn = W;                    % contagious period

% clinical parameters
%--------------------------------------------------------------------------
C.Tic = W;                    % time until symptoms
C.sev = W;                    % P(severe symptoms | symptomatic)
C.Tsy = W;                    % symptomatic period
C.Trd = W;                    % period in CCU
C.fat = W;                    % P(fatality | CCU)
C.sur = W;                    % P(fatality | home)

% testing parameters
%--------------------------------------------------------------------------
C.ttt = U;                    % test, track and trace
C.ont = U;                    % testing latency (months)
C.del = W;                    % test delay (days)
C.tes = V;                    % test selectivity (for infection)
C.sus = U;                    % sustained testing
C.bas = V;                    % baseline testing

% immunity
%--------------------------------------------------------------------------
C.Tim = W;                    % period of immunity
C.Tex = W;                    % period of exemption
C.r   = W;                    % proportion of people not susceptible
C.res = W;                    % proportion with innate immunity
C.tts = U;                    % testing buildup



% log transform
%==========================================================================
P         = spm_vecfun(P,@log);

% field names of random effects
%--------------------------------------------------------------------------
field     = fieldnames(P);
str.field = field(rfx);

return