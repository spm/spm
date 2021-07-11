function [P,C,str] = spm_SARS_priors(nN)
% Generate prior expectation and covariance log parameters
% FORMAT [pE,pC,str,rfx] = spm_SARS_priors(nN)
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
% $Id: spm_SARS_priors.m 8123 2021-07-11 10:28:01Z karl $

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
    [P,C,str] = spm_SARS_priors;
    free  = {'N','Nin','Nou','sde','mem','qua','hos','ccu','res','Tim',...
             'sev','lat','fat','sur','tes','tts','rol','pro'};

    if nN == 1
        
        return
        
    elseif nN == 3
        
        % free parameters for mixture model (age groups)
        %------------------------------------------------------------------
        for i = 1:numel(free)
            P.(free{i}) = kron(ones(nN,1),P.(free{i}));
            C.(free{i}) = kron(ones(nN,1),C.(free{i}));
        end
        
        % Age-specific prior expectations
        %------------------------------------------------------------------
        P.N   = P.N - log(nN);
        P.n   = P.n - log(nN);
        P.rol = log([0.01  (365 + 0)   32;
                     0.02  (365 + 0)   32;
                     0.04  (365 + 0)   32]);
        P.fol = log([0.02  (365 + 128) 32;
                     0.02  (365 + 64)  32;
                     0.02  (365 + 32)  32]);
                 
        P.sev = log([0.00002;
                     0.002;
                     0.03]);
                 
        P.lat = P.sev;
        
        % contact matrices: number of contacts per day
        %------------------------------------------------------------------
        P.Nin = log([2     1     1;
                     1     2     1;
                     1     1     2]/2);
        
        P.Nou = log([32    4     1;
                     4     24    4;
                     1     4     16]);
        
        C.Nin = (exp(P.Nin) > 0)*C.Nin(1);
        C.Nou = (exp(P.Nou) > 0)*C.Nou(1);
        
        return
        
    elseif nN == 4
        
        % free parameters for mixture model (age groups)
        %----------------------------------------------------------------------
        for i = 1:numel(free)
            P.(free{i}) = kron(ones(nN,1),P.(free{i}));
            C.(free{i}) = kron(ones(nN,1),C.(free{i}));
        end
        
        % Age-specific prior expectations
        %----------------------------------------------------------------------
        P.N   = P.N - log(nN);
        P.n   = P.n - log(nN);
        P.rol = log([0.0001 (365 + 0)   32;
                     0.001  (365 + 0)   32;
                     0.01   (365 + 0)   32;
                     0.02   (365 + 0)   32]);
        P.fol = log([0.01   (365 + 512) 32;
                     0.01   (365 + 128) 32;
                     0.02   (365 + 64 ) 32;
                     0.02   (365 + 32 ) 32]);
                 
        C.rol = [0    0    0   ;
                 1/16 1/64 1/64;
                 1/16 1/64 1/64;
                 1/16 1/64 1/64];
        C.fol = [0    0    0   ;
                 1/16 1/64 1/64;
                 1/16 1/64 1/64;
                 1/16 1/64 1/64];
             
        P.sev = log([0.000005;
                     0.0001;
                     0.002;
                     0.04]);
        P.lat = P.sev;
        
        % unvaccinated proportion
        %------------------------------------------------------------------
        P.pro = log([0.64;
                     0.32;
                     0.16;
                     0.08]);
        
        % contact matrices: number of contacts per day
        %------------------------------------------------------------------
        P.Nin = log([2     1     1    1;
                     1     2     1    1;
                     1     1     2    1;
                     1     1     1    2]/2);
        
        P.Nou = log([24    4     2    1;
                     4     32    4    2;
                     2     4     16   4;
                     1     2     4    16]);
        
        C.Nin = (exp(P.Nin) > 0)*C.Nin(1);
        C.Nou = (exp(P.Nou) > 0)*C.Nou(1);
        
        return
    end
    
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
names{7}  = 'time constant of lockdown';
names{8}  = 'time constant of unlockin';
names{9}  = 'viral spreading (days)';
names{10} = 'admission rate (hospital)';
names{11} = 'admission rate (critical)';
names{12} = 'contact rates exponent';

% infection (transmission) parameters
%--------------------------------------------------------------------------
names{13} = 'contacts: home';
names{14} = 'contacts: work';
names{15} = 'transmission strength';
names{16} = 'seasonal transmission';
names{17} = 'infected period   (days)';
names{18} = 'infectious period (days)';
names{19} = 'loss of immunity  (days)';
names{20} = 'resistance (late)';

% clinical parameters
%--------------------------------------------------------------------------
names{21} = 'asymptomatic period (days)';
names{22} = 'symptomatic period (days)';
names{23} = 'critical period (days)';
names{24} = 'P(ARDS|symptoms): early';
names{25} = 'P(ARDS|symptoms): late';
names{26} = 'P(fatality|ARDS): early';
names{27} = 'P(fatality|ARDS): late';

% testing parameters
%--------------------------------------------------------------------------
names{28} = 'FTTI efficacy';
names{29} = 'testing: bias (PCR)';
names{30} = 'testing: bias (LFD)';
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
names{41} = 'vaccination rollout (1st)';
names{42} = 'vaccination rollout (2nd)';

names{43} = 'vaccine efficacy: sterilising';
names{44} = 'vaccine efficacy: pathogenicity';
names{45} = 'vaccine efficacy: transmission';
names{46} = 'vaccine efficacy: fatality';

names{47} = 'LFD confirmation';
names{48} = 'self-isolation (days)';
names{49} = 'loss of T-cell immunity';

names{50} = 'LFD specificity';
names{51} = 'LFD sensitivity';
names{52} = 'relative positivity';
names{53} = 'unvaccinated proportion';



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
% Y(:,26) - Population immunity (total)
% Y(:,27) - Hospital cases
% Y(:,28) - Incidence of Long Covid
% Y(:,29) - population immunity (vaccine)

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
    'Hospital cases'...
    'Incidence of Long Covid'...
    'Vaccine seroconversion'};

str.factors = factors;
str.factor  = factor;
str.names   = names;

% Expectations (either heuristic or taken from the above sources)
%==========================================================================
P.N   = 64;                   % (01) population size (millions)
P.n   = exp(-8);              % (02) proportion of initial cases (cases)
P.r   = 0.1;                  % (03) pre-existing immunity (proportion)
P.o   = 0.1;                  % (04) initial exposed proportion
P.m   = 0.1;                  % (05) relative eflux

% location (exposure) parameters
%--------------------------------------------------------------------------
P.out = 0.4;                  % (06) P(leaving home)
P.sde = 4;                    % (07) time constant of lockdown
P.qua = 64;                   % (08) time constant of unlocking
P.exp = 0.02;                 % (09) viral spreading (rate)
P.hos = 1;                    % (10) admission rate (hospital) [erf]
P.ccu = 0.2;                  % (11) admission rate (CCU)
P.s   = 2;                    % (12) exponent of contact rates

% infection (transmission) parameters
%--------------------------------------------------------------------------
P.Nin = 2;                    % (13) effective number of contacts: home
P.Nou = 24;                   % (14) effective number of contacts: work
P.trn = 0.2;                  % (15) transmission strength (secondary attack rate)
P.trm = 0.04;                 % (16) seasonality
P.Tin = 3;                    % (17) infected period (days)
P.Tcn = 4;                    % (18) infectious period (days)
P.Tim = 128;                  % (19) seropositive immunity (days)
P.res = 0.2;                  % (20) seronegative proportion (late)

% clinical parameters
%--------------------------------------------------------------------------
P.Tic = 4;                    % (21) asymptomatic period (days)
P.Tsy = 8;                    % (22) symptomatic period  (days)
P.Trd = 6;                    % (23) CCU period (days)

P.sev = 0.04;                 % (24) P(ARDS | symptoms): winter
P.lat = 0.04;                 % (25) P(ARDS | symptoms): summer
P.fat = 0.5;                  % (26) P(fatality | ARDS): winter
P.sur = 0.5;                  % (27) P(fatality | ARDS): summer

% testing parameters
%--------------------------------------------------------------------------
P.ttt = 0.036;                % (28) FTTI efficacy
P.tes = [16 8];               % (29) bias (for infection): PCR (Pill. 1 & 2)
P.tts = 1;                    % (30) bias (for infection): LFD
P.del = 3;                    % (31) test delay (days)
P.vac = 512;                  % (32) vaccination time constant (days)
P.fnr = [0.2 0.1];            % (33) false-negative rate  (infected/ious]
P.fpr = [0.002 0.02];         % (34) false-positive rate: (Sus. and Ab +ve)

P.lim = [1   2   2   2]/1000; % (35) testing: capacity (P1, P2, & LFD)
P.rat = [8   24  8   4];      % (36) testing: dispersion
P.ons = [100 200 300 400];    % (37) testing: onset

P.lag = [1 1];                % (38) reporting lag
P.inn = 1;                    % (39) seasonal phase
P.mem = 256;                  % (40) unlocking time constant
P.rol = [0.02 365 32 0.02];   % (41) vaccination rollout (1st)
P.fol = [0.02 450 32 0.02];   % (42) vaccination rollout (2nd)

P.vef = 0.8;                  % (43) vaccine efficacy: infection
P.lnk = 0.1;                  % (44) vaccine efficacy: pathogenicity
P.ves = 0.1;                  % (45) vaccine efficacy: transmission
P.lnf = 0.1;                  % (46) vaccine efficacy: fatality

P.con = 0.2;                  % (47) LFD confirmation
P.iso = 8;                    % (48) self-isolation (days)
P.Tnn = 512;                  % (49) loss of T-cell immunity (days)

P.lnr = 0.5;                  % (50) LFD sensitivity
P.lpr = 0.001;                % (51) LFD specificity
P.rel = .5;                   % (52) relative positivity
P.pro = .08;                  % (53) unvaccinated proportion


% infection fatality (for susceptible population)
%--------------------------------------------------------------------------
% IFR (hospital): P.sev*P.fat*100
% IFR (hospital): exp(Ep.sev)*exp(Ep.fat)*100

% Variances (mildly informative priors, apart from initial cases and size)
%==========================================================================
U     = exp( 2);              % flat priors
V     = exp(-2);              % uninformative priors
W     = exp(-4);              % informative priors
X     = exp(-6);              % informative priors
Z     = exp(-8);              % informative priors

C.N   = U;                    % (01) population size (millions)
C.n   = U;                    % (02) initial cases (cases)
C.r   = W;                    % (03) pre-existing immunity (proportion)
C.o   = W;                    % (04) initial exposed proportion
C.m   = W;                    % (05) relative eflux

% location (exposure) parameters
%--------------------------------------------------------------------------
C.out = X;                    % (06) P(leaving home)
C.sde = W;                    % (07) sensitivity to susceptibility
C.qua = W;                    % (08) sensitivity to seroprevalence
C.exp = W;                    % (09) viral spreading (rate)
C.hos = W;                    % (10) admission rate (hospital)
C.ccu = W;                    % (11) admission rate (CCU)
C.s   = W;                    % (12) decay of social distancing

% infection (transmission) parameters
%--------------------------------------------------------------------------
C.Nin = W;                    % (13) effective number of contacts: home
C.Nou = W;                    % (14) effective number of contacts: work
C.trn = X;                    % (16) transmission strength
C.trm = W;                    % (15) seasonality
C.Tin = Z;                    % (17) infected period (days)
C.Tcn = Z;                    % (18) infectious period (days)
C.Tim = X;                    % (19) seropositive immunity (months)
C.res = X;                    % (20) seronegative immunity (proportion)

% clinical parameters
%--------------------------------------------------------------------------
C.Tic = Z;                    % (21) asymptomatic period (days)
C.Tsy = Z;                    % (22) symptomatic period  (days)
C.Trd = X;                    % (23) CCU period (days)
C.sev = W;                    % (24) P(ARDS | symptoms): winter
C.lat = W;                    % (25) P(ARDS | symptoms): summer
C.fat = X;                    % (26) P(fatality | ARDS): rate
C.sur = X;                    % (27) P(fatality | ARDS): decay

% testing parameters
%--------------------------------------------------------------------------
C.ttt = X;                    % (28) FTTI efficacy
C.tes = V;                    % (29) testing: bias (early)
C.tts = V;                    % (30) testing: bias (late)
C.del = X;                    % (31) test delay (days)
C.vac = X;                    % (32) vaccination time constant (days)
C.fnr = W;                    % (33) false-negative rate
C.fpr = W;                    % (34) false-positive rate

C.lim = V;                    % (35) testing: capacity
C.rat = X;                    % (36) testing: constant (days)
C.ons = U;                    % (37) testing: onset (days)

C.lag = V;                    % (38) reporting lag
C.inn = V;                    % (39) seasonal phase
C.mem = V;                    % (40) unlocking time constant
C.rol = X;                    % (41) vaccination rollout (1st)
C.fol = X;                    % (42) vaccination rollout (2nd)

C.vef = W;                    % (43) vaccine efficacy: infection
C.lnk = W;                    % (44) vaccine efficacy: pathogenicity
C.ves = W;                    % (45) vaccine efficacy: transmission
C.lnf = W;                    % (46) vaccine efficacy: fatality

C.con = V;                    % (47) LFD confirmation
C.iso = W;                    % (48) self-isolation (days)
C.Tnn = X;                    % (49) loss of T-cell immunity (days)

C.lnr = W;                    % (50) LFD sensitivity
C.lpr = W;                    % (51) LFD specificity
C.rel = V;                    % (52) relative positivity
C.pro = V;                    % (53) unvaccinated proportion

% check prior expectations and covariances are consistent
%--------------------------------------------------------------------------
field = fieldnames(P);
for i = 1:numel(field)
    if numel(P.(field{i})) ~= numel(C.(field{i}))
        C.(field{i}) = repmat(C.(field{i}),1,numel(P.(field{i})));
    end
end

% log transform
%==========================================================================
P   = spm_vecfun(P,@log);

% field names of random effects
%--------------------------------------------------------------------------
str.field = fieldnames(P);

return