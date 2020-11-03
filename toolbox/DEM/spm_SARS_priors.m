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
% $Id: spm_SARS_priors.m 8001 2020-11-03 19:05:40Z karl $

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
% https://royalsociety.org/-/media/policy/projects/set-c/set-covid-19-R-estimates.pdf
% https://arxiv.org/abs/2006.01283
%--------------------------------------------------------------------------

% parameter names (where %** denotes fixed effects)
%==========================================================================
names{1}  = 'population size';
names{2}  = 'initial cases';
names{3}  = 'pre-existing immunity';
names{4}  = 'initial exposed';


% location (exposure) parameters
%--------------------------------------------------------------------------
names{5}  = 'P(leaving home)';
names{6}  = 'threshold: lockdown';  
names{7}  = 'threshold: containment';
names{8}  = 'viral spreading';
names{9}  = 'bed availability';
names{10} = 'distancing sensitivity';
names{11} = 'quarantine sensitivity';
names{12} = 'mechanical sensitivity';

% infection (transmission) parameters
%--------------------------------------------------------------------------
names{13} = 'contacts: home';
names{14} = 'contacts: work';
names{15} = 'transmission (max)';
names{16} = 'transmission (min)';
names{17} = 'infected period';
names{18} = 'infectious period';
names{19} = 'loss of immunity';
names{20} = 'seronegative proportion';

% clinical parameters
%--------------------------------------------------------------------------
names{21} = 'incubation period';
names{22} = 'symptomatic period';
names{23} = 'hospitalisation period';
names{24} = 'P(ARDS|symptoms)';
names{25} = 'P(fatality|ARDS)';
names{26} = 'P(ARDS|symptoms): days';
names{27} = 'P(fatality|ARDS): days';

% testing parameters
%--------------------------------------------------------------------------
names{28} = 'FTTI efficacy';
names{29} = 'infection bias';
names{30} = 'test delay (days)';
names{31} = 'seropositive-dependent';
names{32} = 'symptom-dependent';
names{33} = 'false-negative rate';
names{34} = 'false-positive rate';

names{35} = 'testing: bias';
names{36} = 'testing: capacity';
names{37} = 'testing: constant';
names{38} = 'testing: onset';

names{39} = 'relative eflux';

% latent or hidden factors
%--------------------------------------------------------------------------
factors   = {'Location','Infection','Symptoms','Testing'};

factor{1} = {'lo-risk','hi-risk','ccu','no-risk','isolated'};
factor{2} = {'susceptible','infected','infectious','AB +ve','AB -ve'};
factor{3} = {'none','symptoms','severe','deceased'};
factor{4} = {'untested','waiting','PCR +ve','PCR -ve'};

% labels or strings for plotting
%--------------------------------------------------------------------------
% Y(:,1)  - daily deaths
% Y(:,2)  - daily tests
% Y(:,3)  - CCU occupancy
% Y(:,4)  - reproduction ratio (R)
% Y(:,5)  - seropositive immunity (%)
% Y(:,6)  - PCR testing rate
% Y(:,7)  - contagion risk (%)
% Y(:,8)  - prevalence (%)
% Y(:,9)  - daily contacts
% Y(:,10) - daily incidence
% Y(:,11) - number infected  
% Y(:,12) - number symptomatic
% Y(:,13) - mobility (%)
% Y(:,14) - work (%)
% Y(:,15) - hospital admissions
% Y(:,16) - fatality rates
% Y(:,17) - home (%)
str.outcome = {'Daily deaths',...
               'Daily tests',...
               'CCU occupancy',...
               'Reproduction ratio',...
               'Seropositive immunity',...
               'PCR testing rate',...
               'Contagion risk (%)',...
               'Prevalence {%}',...
               'Daily contacts',...
               'daily incidence',...
               'Number infected'...
               'Number symptomatic'...
               'Mobility (%)'...
               'Workplace (%)'...
               'Admissions'...
               'Death rates'...
               'Residential (%)'};
str.factors = factors;
str.factor  = factor;
str.names   = names;

% cut and paste to see the effects of changing different prior expectations
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if false
    pE    = spm_SARS_priors; M.T = 12*32;
    [Y,X] = spm_SARS_gen(pE,M,1:3);
    u     = exp(pE.cap + pE.N)*1e6/2;
    spm_SARS_plot(Y,X,[],u)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% Expectations (either heuristic or taken from the above sources)
%==========================================================================
P.N   = 64;                   % (01) population size (millions)
P.n   = 1;                    % (02) initial cases (cases)
P.r   = 0.2;                  % (03) pre-existing immunity (proportion)
P.o   = 0.2;                  % (04) initial exposed proportion

% location (exposure) parameters
%--------------------------------------------------------------------------
P.out = 0.5;                  % (05) P(leaving home)
P.sde = 0.03;                 % (06) lockdown threshold
P.qua = 0.004;                % (07) Quarantine threshold
P.exp = 0.0008;               % (08) P(leaving area)
P.cap = 16/100000;            % (09) bed availability (per capita)
P.s   = 2;                    % (10) distancing sensitivity
P.u   = 6;                    % (11) quarantine sensitivity
P.c   = 1;                    % (12) mechanical sensitivity

% infection (transmission) parameters
%--------------------------------------------------------------------------
P.Nin = 2;                    % (13) effective number of contacts: home
P.Nou = 64;                   % (14) effective number of contacts: work
P.Nru = 0.5;                  % (15) transmission strength (max)
P.trn = 0.5;                  % (16) transmission strength (min)
P.Tin = 5;                    % (17) infected period (days)
P.Tcn = 5;                    % (18) infectious period (days)
P.Tim = 256;                  % (19) seropositive immunity (days)
P.res = 0.5;                  % (20) seronegative proportion

% clinical parameters
%--------------------------------------------------------------------------
P.Tic = 4;                    % (21) incubation period (days)
P.Tsy = 6;                    % (22) symptomatic period (days)
P.Trd = 10;                   % (23) CCU period (days)

P.sev = 0.3/100;              % (24) P(ARDS | symptoms)
P.fat = 0.5;                  % (25) P(fatality | ARDS/symptoms)
P.lat = 0.2/100;              % (26) P(ARDS | symptoms): late
P.sur = 0.5;                  % (27) P(fatality | ARDS): late

% testing parameters
%--------------------------------------------------------------------------
P.ttt = 0.036;                % (28) FTTI efficacy
P.tes = 2;                    % (29) test selectivity (for infection)
P.del = 4;                    % (30) test delay (days)
P.sus = 0.01;                 % (31) seropositive-dependent
P.ont = 0.01;                 % (32) symptom-dependent
P.fnr = 0.2;                  % (33) false-negative rate
P.fpr = 0.002;                % (34) false-positive rate

P.tts = 256;                  % (35) testing: bias
P.lim = [0.02 0.002];         % (36) testing: capacity
P.rat = [64    4];            % (37) testing: constant
P.ons = [256 100];            % (38) testing: onset

P.m   = 4;                    % (39) relative eflux



% cut and paste to see the effects of changing different prior expectations
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if false
    pE    = spm_SARS_priors; M.T = 12*32;
    [Y,X] = spm_SARS_gen(pE,M,1:3);
    u     = exp(pE.cap + pE.N)*1e6/2;
    spm_SARS_plot(Y,X,[],u)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% infection fatality (for susceptible population)
%--------------------------------------------------------------------------
% IFR (hospital): P.sev*P.fat*100
% IFR (hospital): exp(Ep.sev)*exp(Ep.fat)*100

% Variances (mildly informative priors, apart from initial cases and size)
%==========================================================================
U     = 4;                    % flat priors
V     = 1/16;                 % uninformative priors
W     = 1/64;                 % informative priors
X     = 1/1024;               % informative priors

C.N   = U;                    % (01) population size (millions)
C.n   = U;                    % (02) initial cases (cases)
C.r   = V;                    % (03) pre-existing immunity (proportion)
C.o   = V;                    % (04) initial exposed proportion

% location (exposure) parameters
%--------------------------------------------------------------------------
C.out = X;                    % (05) P(leaving home)
C.sde = V;                    % (06) lockdown threshold
C.qua = V;                    % (07) Quarantine threshold
C.exp = V;                    % (08) contacts (viral spreading)
C.cap = V;                    % (09) bed availability (per capita)
C.s   = W;                    % (10) distancing sensitivity
C.u   = W;                    % (11) quarantine sensitivity 
C.c   = W;                    % (12) mechanical sensitivity

% infection (transmission) parameters
%--------------------------------------------------------------------------
C.Nin = V;                    % (13) effective number of contacts: home
C.Nou = V;                    % (14) effective number of contacts: work
C.Nru = W;                    % (15) transmission strength (max)
C.trn = W;                    % (16) transmission strength (min)
C.Tin = X;                    % (17) infected period (days)
C.Tcn = X;                    % (18) infectious period (days)
C.Tim = W;                    % (19) seropositive immunity (months)
C.res = W;                    % (20) seronegative immunity (proportion)

% clinical parameters
%--------------------------------------------------------------------------
C.Tic = X;                    % (21) incubation period  (days)
C.Tsy = X;                    % (22) symptomatic period (days)
C.Trd = X;                    % (23) CCU period (days)
C.sev = X;                    % (24) P(ARDS | symptoms)
C.fat = X;                    % (25) P(fatality | ARDS)
C.lat = U;                    % (26) P(ARDS | symptoms): decay
C.sur = U;                    % (27) P(fatality | ARDS): decay

% testing parameters
%--------------------------------------------------------------------------
C.ttt = X;                    % (28) FTTI efficacy
C.tes = U;                    % (29) test selectivity (for infection)
C.del = X;                    % (30) test delay (days)
C.sus = U;                    % (31) seropositive-dependent
C.ont = U;                    % (32) symptom-dependent
C.fnr = X;                    % (33) false-negative rate
C.fpr = X;                    % (34) false-positive rate

C.tts = U;                    % (25) testing: bias
C.lim = U;                    % (36) testing: capacity
C.rat = U;                    % (37) testing: constant (days)
C.ons = U;                    % (38) testing: onset (days)

C.m   = W;                    % (39) relative eflux


% implicit prior confidence bounds
%--------------------------------------------------------------------------
s     = spm_invNcdf(1 - 0.01);
field = fieldnames(P);
for i = 1:numel(field)
    b = log(P.(field{i})) + [-s s]*sqrt(C.(field{i}));
    R.(field{i}) = exp(b);
    
    % check prior expectations and covariances are consistent
    %---------------------------------------------------------------------
    if numel(P.(field{i})) ~= numel(C.(field{i}))
        C.(field{i}) = repmat(C.(field{i}),1,numel(P.(field{i})));
    end
end

% log transform
%==========================================================================
P   = spm_vecfun(P,@log);

% cut and paste to see the effects of changing different prior expectations
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if false
    pE    = spm_SARS_priors; M.T = 12*32;
    [Y,X] = spm_SARS_gen(pE,M,1:3);
    u     = exp(pE.cap + pE.N)*1e6/2;
    spm_SARS_plot(Y,X,[],u)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




% field names of random effects
%--------------------------------------------------------------------------
str.field = fieldnames(P);

return