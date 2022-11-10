function [P,C,str] = spm_SARS_priors_B(nN)
% Generate prior expectation and covariance log parameters (bound version)
% FORMAT [pE,pC,str] = spm_SARS_priors_B(nN)
%
% nN          - number of age groups
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
% https://royalsociety.org/-/media/policy/projects/set-c/set-covid-19-R-estimates.pdf
% https://arxiv.org/abs/2006.01283
%--------------------------------------------------------------------------

% priors for multiple age groups
%==========================================================================
if nargin
    
    % priors for single group
    %----------------------------------------------------------------------
    [P,C,str] = spm_SARS_priors_B;
    
    if nN == 1, return, end
    
    % free parameters for mixture model (age groups)
    %----------------------------------------------------------------------
    free  = {'N','Nin','Nou','hos','ccu','res','sev','lat','fat','sur','tes','tts','rol'};
    for i = 1:numel(free)
        P.(free{i}) = kron(ones(nN,1),P.(free{i}));
        C.(free{i}) = kron(ones(nN,1),C.(free{i}));
    end
    
    % Age-specific prior expectations
    %----------------------------------------------------------------------
    P.N   = P.N - log(nN);
    P.n   = P.n - 8;
    P.rol = log([0.01  (365 + 256) 128 0.01;
                 0.02  (365 + 64 ) 64  0.02;
                 0.04  (365 + 0  ) 16  0.04]);
    C.rol = [1/16 1/64 1/64 1/16;
             1/16 1/64 1/64 1/16;
             1/16 1/64 1/64 1/16];
    P.sev = log([0.00002;
                 0.002;
                 0.02]);
    P.lat = P.sev;
    
    % contact matrices: number of contacts per day
    %----------------------------------------------------------------------
    P.Nin = log([2     1     1;
                 1     2     1;
                 1     1     2]/2);
    
    P.Nou = log([32    4     4;
                 4     24    4;
                 4     4     16]);
             
    C.Nin = (exp(P.Nin) > 0)*C.Nin(1);
    C.Nou = (exp(P.Nou) > 0)*C.Nou(1);
    
    return
end

% parameter names
%==========================================================================
names{1}  = 'population size (M)';
names{2}  = 'initial cases';
names{3}  = 'pre-existing immunity';
names{4}  = 'initially exposed';
names{5}  = 'relative eflux';

% location (exposure) parameters
%--------------------------------------------------------------------------
names{6}  = 'P(leaving home)';
names{7}  = 'sensitivity to susceptible';
names{8}  = 'sensitivity to prevalence';
names{9}  = 'viral spreading (days)';
names{10} = 'admission rate (hospital)';
names{11} = 'admission rate (CCU)';
names{12} = 'decay of social distancing';

% infection (transmission) parameters
%--------------------------------------------------------------------------
names{13} = 'contacts: home';
names{14} = 'contacts: work';
names{15} = 'transmission strength';
names{16} = 'seasonality';
names{17} = 'infected period   (days)';
names{18} = 'infectious period (days)';
names{19} = 'loss of immunity  (days)';
names{20} = 'resistance (late)';

% clinical parameters
%--------------------------------------------------------------------------
names{21} = 'asymptomatic period (days)';
names{22} = 'symptomatic period (days)';
names{23} = 'ARDS period (days)';
names{24} = 'P(ARDS|symptoms): winter';
names{25} = 'P(ARDS|symptoms): summer';
names{26} = 'P(fatality|ARDS): rate';
names{27} = 'P(fatality|ARDS): days';

% testing parameters
%--------------------------------------------------------------------------
names{28} = 'FTTI efficacy';
names{29} = 'testing: bias (early)';
names{30} = 'testing: bias (late)';
names{31} = 'test delay (days)';
names{32} = 'vaccine constant (days)';
names{33} = 'false-negative rate';
names{34} = 'false-positive rate';

names{35} = 'testing: capacity';
names{36} = 'testing: constant';
names{37} = 'testing: onset';

names{38} = 'reporting lag';
names{39} = 'seasonal phase';
names{40} = 'unlocking time constant';
names{41} = 'vaccination rollout';
names{42} = 'loss of T-cell immunity';
names{43} = 'vaccine efficacy: infection';
names{44} = 'vaccine efficacy: pathogenicity';
names{45} = 'vaccine efficacy: transmission';


% latent or hidden factors
%--------------------------------------------------------------------------
factors   = {'Location','Infection','Symptoms','Testing'};

factor{1} = {'lo-risk','hi-risk','ICU','no-risk','isolated','hospital'};
factor{2} = {'susceptible','infected','infectious','Ab +ve','Ab -ve','Vac +ve'};
factor{3} = {'none','symptoms','severe','deceased'};
factor{4} = {'untested','waiting','PCR +ve','PCR -ve','LFD +ve','LFD -ve'};
factor{5} = {' ',' '};

% labels or strings for plotting
%--------------------------------------------------------------------------
% Y(:,1)  - Daily deaths (28 days)
% Y(:,2)  - Daily confirmed cases
% Y(:,3)  - Mechanical ventilation
% Y(:,4)  - Reproduction ratio (R)
% Y(:,5)  - Seroprevalence {%}
% Y(:,6)  - PCR testing rate
% Y(:,7)  - Contagion risk (%)
% Y(:,8)  - Contagious {%}
% Y(:,9)  - Daily contacts
% Y(:,10) - Daily incidence (%)
% Y(:,11) - Prevalence {%}
% Y(:,12) - Number symptomatic
% Y(:,13) - Mobility (%)
% Y(:,14) - Workplace (%)
% Y(:,15) - Certified deaths
% Y(:,16) - Hospital admissions
% Y(:,17) - Hospital deaths
% Y(:,18) - Non-hospital deaths
% Y(:,19) - Daily incidence (per hundred thousand)
% Y(:,20) - Weekly confirmed cases (per hundred thousand)
% Y(:,21) - Infection fatality ratio (%)
% Y(:,22) - Cumulative first dose
% Y(:,23) - PCR case positivity (%)
% Y(:,24) - Lateral flow tests
% Y(:,25) - Cumulative attack rate
% Y(:,26) - Population immunity
% Y(:,27) - Hospital cases

str.outcome = {'Daily deaths (28 days)',...
    'Daily confirmed cases',...
    'Mechanical ventilation',...
    'Reproduction ratio',...
    'Seroprevalence {%}',...
    'PCR testing rate',...
    'Contagion risk (%)',...
    'Contagious {%}',...
    'Daily contacts',...
    'Daily incidence (%)',...
    'Prevalence (%)'...
    'Number symptomatic'...
    'Mobility (%)'...
    'Workplace (%)'...
    'Certified deaths',...
    'Hospital admissions'...
    'Hospital deaths',...
    'Non-hospital deaths'...
    'Daily incidence (per 100,000)',...
    'Weekly confirmed cases (per 100,000)',...
    'IFR (%)',...
    'Cumulative first dose',...
    'PCR positivity (%)',...
    'Lateral flow tests',...
    'Attack rate (%)',...
    'Herd immunity (%)'...
    'Hospital cases'};

str.factors = factors;
str.factor  = factor;
str.names   = names;

% Expectations (either heuristic or taken from the above sources)
%==========================================================================
P.N   = [4; 1024];             % (01) population size (millions)
P.n   = [exp(-16); 1];         % (02) proportion of initial cases (cases)
P.r   = [4; 24]/100;           % (03) pre-existing immunity (proportion)
P.o   = [4; 24]/100;           % (04) initial exposed proportion
P.m   = [4; 24]/100;           % (05) relative eflux

% location (exposure) parameters
%--------------------------------------------------------------------------
P.out = [0.3; 0.5];            % (06) P(leaving home)
P.sde = [4; 32];               % (07) time constant of lockdown
P.qua = [32; 128];             % (08) time constant of unlocking
P.exp = [1/64; 1/8];           % (09) viral spreading (rate)
P.hos = [0.5; 2];              % (10) admission rate (hospital)
P.ccu = [0.1; 0.4];            % (11) admission rate (CCU)
P.s   = [1/4; 4];              % (12) exponent of contact rates

% infection (transmission) parameters
%--------------------------------------------------------------------------
P.Nin = [1/2; 4];              % (13) effective number of contacts: home
P.Nou = [8; 64];               % (14) effective number of contacts: work
P.trn = [8; 32]/100;           % (15) transmission strength (secondary attack rate)
P.trm = [4; 16]/100;           % (16) seasonality
P.Tin = [3; 5];                % (17) infected period (days)
P.Tcn = [3; 5];                % (18) infectious period (days)
P.Tim = [128; 256];            % (19) seropositive immunity (days)
P.res = [8; 32]/100;           % (20) seronegative proportion (late)

% clinical parameters
%--------------------------------------------------------------------------
P.Tic = [3 3; 5 5];            % (21) asymptomatic period (days)
P.Tsy = [6; 12];               % (22) symptomatic period  (days)
P.Trd = [6; 8];                % (23) CCU period (days)

P.sev = [1/8; 2]/100;          % (24) P(ARDS | symptoms): winter
P.lat = [1/8; 2]/100;          % (25) P(ARDS | symptoms): summer
P.fat = [0.2; 0.8];            % (26) P(fatality | ARDS): winter
P.sur = [0.2; 0.8];            % (27) P(fatality | ARDS): summper

% testing parameters
%--------------------------------------------------------------------------
P.ttt = [1; 4]/100;            % (28) FTTI efficacy
P.tes = [8  4;
         32 16];               % (29) bias (for infection): PCR (Pill. 1 & 2)
P.tts = [1/2; 2];              % (30) bias (for infection): LFD
P.del = [2 4];                 % (31) test delay (days)
P.vac = [256; 1024];           % (32) vaccination time constant (days)
P.fnr = [0.1 0.01; 
         0.4 0.2];             % (33) false-negative rate  (infected/ious]
P.fpr = [0.1 1;
         0.4 4]/100;           % (34) false-positive rate: (Sus. and Ab +ve)

P.lim = [1/8 1   2   2;
         2   4   4   8]/1000;  % (35) testing: capacity (P1, P2, & LFD)
P.rat = [6   16  6   2;
         12  32  12  8];       % (36) testing: dispersion
P.ons = [50  100 200 300;
         100 200 300 500];     % (37) testing: onset

P.lag = [1  1;
         16 16];               % (38) reporting lag
P.inn = [1; 32];               % (39) seasonal phase
P.mem = [128; 512];            % (40) unlocking time constant
P.rol = [0.01 365 32 0.01;
         0.08 512 64 0.02];    % (41) vaccination rollout
P.Tnn = [512; 1024];           % (42) loss of T-cell immunity (days)
P.vef = [0.1; 0.3];            % (43) vaccine efficacy: infection
P.lnk = [0.1; 0.3];            % (44) vaccine efficacy: pathogenicity
P.ves = [0.1; 0.3];            % (45) vaccine efficacy: transmission

% Variances and log transform mean
%--------------------------------------------------------------------------
field = fieldnames(P);
for i = 1:numel(field)
    C.(field{i}) = diff(log(P.(field{i})))/32;
    P.(field{i}) = mean(log(P.(field{i})));
end

% field names of random effects
%--------------------------------------------------------------------------
str.field = fieldnames(P);

return