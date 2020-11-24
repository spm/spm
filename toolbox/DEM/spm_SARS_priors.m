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
% $Id: spm_SARS_priors.m 8017 2020-11-24 21:48:26Z karl $

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
names{1}  = 'population size (M)';
names{2}  = 'initial cases';
names{3}  = 'pre-existing immunity';
names{4}  = 'initially exposed';
names{5}  = 'relative eflux';

% location (exposure) parameters
%--------------------------------------------------------------------------
names{6}  = 'P(leaving home)';
names{7}  = 'lockdown threshold';  
names{8}  = 'seropositive contribution';
names{9}  = 'viral spreading (days)';
names{10} = 'admission rate (hospital)';
names{11} = 'admission rate (CCU)';
names{12} = 'distancing sensitivity';

% infection (transmission) parameters
%--------------------------------------------------------------------------
names{13} = 'contacts: home';
names{14} = 'contacts: work';
names{15} = 'transmission (early)';
names{16} = 'transmission (late)';
names{17} = 'infected period   (days)';
names{18} = 'infectious period (days)';
names{19} = 'loss of immunity  (days)';
names{20} = 'seronegative proportion';

% clinical parameters
%--------------------------------------------------------------------------
names{21} = 'incubation period  (days)';
names{22} = 'symptomatic period (days)';
names{23} = 'ARDS period (days)';
names{24} = 'P(ARDS|symptoms): early';
names{25} = 'P(ARDS|symptoms): late';
names{26} = 'P(fatality|ARDS): early';
names{27} = 'P(fatality|ARDS): late';

% testing parameters
%--------------------------------------------------------------------------
names{28} = 'FTTI efficacy';
names{29} = 'testing: bias (early)';
names{30} = 'testing: bias (late)';
names{31} = 'test delay (days)';
names{32} = 'symptom demand';
names{33} = 'false-negative rate';
names{34} = 'false-positive rate';

names{35} = 'testing: capacity';
names{36} = 'testing: constant';
names{37} = 'testing: onset';

names{38} = 'ventilation (days)';


% latent or hidden factors
%--------------------------------------------------------------------------
factors   = {'Location','Infection','Symptoms','Testing'};

factor{1} = {'lo-risk','hi-risk','ccu','no-risk','isolated','hospital'};
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
% Y(:,9)  - new contacts per day
% Y(:,10) - daily incidence
% Y(:,11) - number infected  
% Y(:,12) - number symptomatic
% Y(:,13) - mobility (%)
% Y(:,14) - work (%)
% Y(:,15) - certified deaths/day
% Y(:,16) - hospitalisation
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
               'Certified deaths',...
               'Admissions'};
           
str.factors = factors;
str.factor  = factor;
str.names   = names;

% cut and paste to see the effects of changing different prior expectations
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if false
    U     = [1 2 3 16];
    pE    = spm_SARS_priors; M.T = 10*32;
    [Y,X] = spm_SARS_gen(pE,M,U);
    spm_SARS_plot(Y,X,[],[],U)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% Expectations (either heuristic or taken from the above sources)
%==========================================================================
P.N   = 64;                   % (01) population size (millions)
P.n   = 1;                    % (02) initial cases (cases)
P.r   = 0.4;                  % (03) pre-existing immunity (proportion)
P.o   = 0.02;                 % (04) initial exposed proportion
P.m   = 0.3;                  % (05) relative eflux

% location (exposure) parameters
%--------------------------------------------------------------------------
P.out = 0.5;                  % (06) P(leaving home)
P.sde = 0.03;                 % (07) lockdown threshold
P.qua = 0.6;                  % (08) seropositive contribution
P.exp = 0.0065;               % (09) viral spreading (days)
P.hos = 0.72;                 % (10) admission rate (hospital)
P.ccu = 0.01;                 % (11) admission rate (CCU)
P.s   = 2.5;                  % (12) distancing sensitivity

% infection (transmission) parameters
%--------------------------------------------------------------------------
P.Nin = 1.5;                  % (13) effective number of contacts: home
P.Nou = 48;                   % (14) effective number of contacts: work
P.trn = 0.4;                  % (15) transmission strength (early)
P.trm = 0.2;                  % (16) transmission strength (late)
P.Tin = 6;                    % (17) infected period (days)
P.Tcn = 5;                    % (18) infectious period (days)
P.Tim = 200;                  % (19) seropositive immunity (days)
P.res = 0.12;                 % (20) seronegative proportion

% clinical parameters
%--------------------------------------------------------------------------
P.Tic = 6;                    % (21) incubation period (days)
P.Tsy = 8;                    % (22) symptomatic period (days)
P.Trd = 12;                   % (23) CCU period (days)

P.sev = 0.5/100;              % (24) P(ARDS | symptoms): early
P.lat = 0.4/100;              % (25) P(ARDS | symptoms): late
P.fat = 0.5;                  % (26) P(fatality | ARDS): early
P.sur = 0.4;                  % (27) P(fatality | ARDS): late

% testing parameters
%--------------------------------------------------------------------------
P.ttt = 0.036;                % (28) FTTI efficacy
P.tes = 6;                    % (29) bias (for infection): early
P.tts = 8;                    % (30) bias (for infection): late
P.del = 4;                    % (31) test delay (days)
P.ont = 0.0002;               % (32) symptom-dependent
P.fnr = 0.2;                  % (33) false-negative rate
P.fpr = 0.002;                % (34) false-positive rate

P.lim = [0.004 0.001];        % (35) testing: capacity
P.rat = [48    8];            % (36) testing: dispersion
P.ons = [200 100];            % (37) testing: onset

P.tcu = 256;                  % (38) ventilation (days)



% cut and paste to see the effects of changing different prior expectations
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if false
    U     = [1 2 3 16];
    pE    = spm_SARS_priors; M.T = 10*32;
    [Y,X] = spm_SARS_gen(pE,M,U);
    spm_SARS_plot(Y,X,[],[],U)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% infection fatality (for susceptible population)
%--------------------------------------------------------------------------
% IFR (hospital): P.sev*P.fat*100
% IFR (hospital): exp(Ep.sev)*exp(Ep.fat)*100

% Variances (mildly informative priors, apart from initial cases and size)
%==========================================================================
U     = exp( 4);              % flat priors
V     = exp(-2);              % uninformative priors
W     = exp(-4);              % informative priors
X     = exp(-6);              % informative priors

C.N   = U;                    % (01) population size (millions)
C.n   = U;                    % (02) initial cases (cases)
C.r   = W;                    % (03) pre-existing immunity (proportion)
C.o   = U;                    % (04) initial exposed proportion
C.m   = V;                    % (05) relative eflux

% location (exposure) parameters
%--------------------------------------------------------------------------
C.out = X;                    % (06) P(leaving home)
C.sde = V;                    % (07) lockdown threshold
C.qua = V;                    % (08) seropositive contribution
C.exp = V;                    % (09) viral spreading (days)
C.hos = X;                    % (10) admission rate (hospital)
C.ccu = W;                    % (11) admission rate (CCU)
C.s   = W;                    % (12) distancing sensitivity

% infection (transmission) parameters
%--------------------------------------------------------------------------
C.Nin = W;                    % (13) effective number of contacts: home
C.Nou = W;                    % (14) effective number of contacts: work
C.trn = W;                    % (16) transmission strength (early)
C.trm = W;                    % (15) transmission strength (late)
C.Tin = X;                    % (17) infected period (days)
C.Tcn = X;                    % (18) infectious period (days)
C.Tim = W;                    % (19) seropositive immunity (months)
C.res = V;                    % (20) seronegative immunity (proportion)

% clinical parameters
%--------------------------------------------------------------------------
C.Tic = X;                    % (21) incubation period  (days)
C.Tsy = X;                    % (22) symptomatic period (days)
C.Trd = X;                    % (23) CCU period (days)
C.sev = X;                    % (24) P(ARDS | symptoms): early
C.lat = X;                    % (25) P(ARDS | symptoms): late
C.fat = X;                    % (26) P(fatality | ARDS): early
C.sur = X;                    % (27) P(fatality | ARDS): late

% testing parameters
%--------------------------------------------------------------------------
C.ttt = X;                    % (28) FTTI efficacy
C.tes = V;                    % (29) testing: bias (early)
C.tts = V;                    % (30) testing: bias (late)
C.del = X;                    % (31) test delay (days)
C.ont = U;                    % (32) symptom-dependent
C.fnr = X;                    % (33) false-negative rate
C.fpr = X;                    % (34) false-positive rate

C.lim = V;                    % (35) testing: capacity
C.rat = U;                    % (36) testing: constant (days)
C.ons = U;                    % (37) testing: onset (days)

C.tcu = W;                    % (38) ventilation (days)


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