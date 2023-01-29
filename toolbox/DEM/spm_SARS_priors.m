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
% https://pubmed.ncbi.nlm.nih.gov/33948610/
%--------------------------------------------------------------------------

% priors for multiple age groups
%==========================================================================
if nargin
    
    % priors for single group              % age-dependent variability in:
    %----------------------------------------------------------------------
    [P,C,str] = spm_SARS_priors;
    free  = {'N'  ,'trm',...               % population and seasonality
             'Tim','Tnn',...               % immunity
             'Nin','Nou',...               % number of contacts
             'rut','out',...               % social distancing
             'hos','ccu',...               % hospitalisation
             'sev','fat',...               % vulnerability
             'tts','lnk',...               % testing and efficacy
             'rol','fol'};                 % vaccination 

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
        P.rol = log([0.004  64;
                     0.04   32;
                     0.04   16]);
        P.fol = P.rol;
                 
        P.sev = log([0.0010;
                     0.0100;
                     0.1000]);

        
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
        
        % Age-specific prior expectations (amplitude and latency)
        %------------------------------------------------------------------
        P.N   = P.N - log(nN);
        P.n   = P.n - log(nN);
        P.rol = log([0.0001 256;
                     0.004  128;
                     0.04   64;
                     0.08   32]);
                 
        P.fol = log([0.0001 64;
                     0.004  48;
                     0.04   32;
                     0.08   16]);

        C.rol = [1/16 1/64;
                 1/16 1/64;
                 1/16 1/64;
                 1/16 1/64];
             
        C.fol = C.rol;
             
        % probability of hospitalisation when seriously ill
        %------------------------------------------------------------------
        P.hos = log([2
                     2
                     1
                     0.4]);
                 
        % probability of transfer to CCU
        %------------------------------------------------------------------
        P.ccu = log([0.2
                     0.04
                     0.02
                     0.005]);
                 
        % morbidity
        %------------------------------------------------------------------
        P.sev = log([0.002
                     0.020
                     0.020
                     0.200]);
        % mortality
        %------------------------------------------------------------------
        P.fat = log([0.000001;
                     0.0001;
                     0.1;
                     0.4]);
        
        % vaccination link (pathogenicity)
        %------------------------------------------------------------------
        P.lnk = log([0.24;
                     0.24;
                     0.48;
                     0.48]);
                 
        C.lnk = [0;
                 exp(-8);
                 exp(-8);
                 exp(-8)];
                 
        
        % contact matrices: number of contacts per day
        %------------------------------------------------------------------
        P.Nin = log([2     1     1    1;
                     1     2     1    1;
                     1     1     2    1;
                     1     1     1    2]/2);
        
        P.Nou = log([16    4     2    1;
                     4     24    4    2;
                     2     4     16   4;
                     1     2     4    8]);
        
        C.Nin = (exp(P.Nin) > 0)*C.Nin(1);
        C.Nou = (exp(P.Nou) > 0)*C.Nou(1);
        
        % sensitivity to prevalence
        %------------------------------------------------------------------
        P.sde = log(ones(1,nN)*4);
        C.sde = ones(1,nN)*C.sde(1);
        
%         % contact matrices: number of contacts per day
%         % https://www.researchgate.net/figure/Social-contact-matrices-Values-and-colours-show-the-mean-number-of-contacts-per-day_fig4_221698133
%         %------------------------------------------------------------------
%         P.Nin = log([1.9    1.4   3.7   0.1;
%                      0.5    4.7   3.3   0.4;
%                      0.4    0.8   3.3   0.5;
%                      0.1    0.2   1.6   0.7]'/2);
%         
%         P.Nou = log([6.5    1.8   7.0   0.1;
%                      1.0   11.6  11.7   0.6;
%                      1.3    2.3  15.0   1.1;
%                      0.1    0.4   5.7   1.4]');
%         
%         C.Nin = (exp(P.Nin) > 0)*exp(-6);
%         C.Nou = (exp(P.Nou) > 0)*exp(-6);
        
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
names{12} = 'infectivity changes';

% infection (transmission) parameters
%--------------------------------------------------------------------------
names{13} = 'contacts: home';
names{14} = 'contacts: work';
names{15} = 'transmission strength';
names{16} = 'seasonal transmission';
names{17} = 'infected period   (days)';
names{18} = 'infectious period (days)';
names{19} = 'loss of natural immunity (days)';
names{20} = 'resistance';

% clinical parameters
%--------------------------------------------------------------------------
names{21} = 'asymptomatic period (days)';
names{22} = 'symptomatic period (days)';
names{23} = 'critical period (days)';
names{24} = 'P(ARDS|symptoms): initial';
names{25} = 'P(ARDS|symptoms): change';
names{26} = 'P(fatality|ARDS): initial';
names{27} = 'P(fatality|ARDS): change';

% testing parameters
%--------------------------------------------------------------------------
names{28} = 'FTTI efficacy';
names{29} = 'testing: bias (PCR)';
names{30} = 'testing: bias (LFD)';
names{31} = 'test delay (days)';
names{32} = 'vaccine seroconversion (days)';
names{33} = 'false-negative rate';
names{34} = 'false-positive rate';

names{35} = 'testing: capacity';
names{36} = 'testing: constant';
names{37} = 'testing: onset';

names{38} = 'reporting lag';
names{39} = 'seasonal phase';
names{40} = 'vaccination rollout (days)';
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
names{52} = 'PCR testing of fatalities';
names{53} = 'contact rate decay (days)';
names{54} = 'survival risk in care homes';
names{55} = 'changes in transfer to CCU';  
names{56} = 'transmissibility parameters';
names{57} = 'doses per seroconversion';
names{58} = 'age-related testing';
names{59} = 'self-isolation';
names{60} = 'Sensitivity to contact rate';


% latent or hidden factors
%--------------------------------------------------------------------------
factors   = {'Location','Infection','Symptoms','Testing'};

factor{1} = {'lo-risk','hi-risk','CCU','safe','isolated','hospital'};
factor{2} = {'susceptible','exposed','infectious','Ab +ve','Ab -ve','vaccinated','exposed (Ab)','infectious (Ab)'};
factor{3} = {'asymptomatic','symptomatic','severe','deceased'};
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
% Y(:,7)  - Risk of infection (%)
% Y(:,8)  - Prevalence (true) {%}
% Y(:,9)  - Daily contacts
% Y(:,10) - Daily incidence (%)
% Y(:,11) - Prevalence (positivity){%}
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
% Y(:,29) - Vaccine antibodies
% Y(:,30) - Cumulative admissions
% Y(:,31) - Vaccine effectiveness (prevalence)
% Y(:,32) - Gross domestic product
% Y(:,33) - Doubling time (days)
% Y(:,34) - Incidence of new cases (total)
% Y(:,35) - Serial interval
% Y(:,36) - Vaccination rate

str.outcome = {'Daily deaths (28 days)',...
    'Daily confirmed cases',...
    'Mechanical ventilation',...
    'Reproduction ratio',...
    'Seroprevalence {%}',...
    'PCR testing rate',...
    'Risk of infection (%)',...
    'Prevalence (true) {%}',...
    'Daily contacts',...
    'Daily incidence (%)',...
    'Prevalence (positivity) (%)'...
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
    'Population immunity (%)'...
    'Hospital cases'...
    'Incidence of Long Covid'...
    'Vaccine antibodies (%)'...
    'Cumulative admissions'...
    'Vaccine effectiveness (prevalence) '...
    'Gross domestic product'...
    'Doubling time (days)'...
    'Daily incidence'...
    'Serial interval'...
    'Vaccination rate'};

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
P.sde = 16;                   % (07) time constant of lockdown
P.qua = 128;                  % (08) time constant of unlocking
P.exp = 0.02;                 % (09) viral spreading (rate)
P.hos = 1;                    % (10) admission rate (hospital) [erf]
P.ccu = 0.1;                  % (11) admission rate (CCU)
P.s   = [1 1 1 1 1 1];        % (12) infectivity exponents

% infection (transmission) parameters
%--------------------------------------------------------------------------
P.Nin = 2;                    % (13) effective number of contacts: home
P.Nou = 24;                   % (14) effective number of contacts: work
P.trn = 0.2;                  % (15) transmission strength (secondary attack rate)
P.trm = 0.04;                 % (16) seasonality
P.Tin = 3;                    % (17) infected period (days)
P.Tcn = 4;                    % (18) infectious period (days)
P.Tim = 128;                  % (19) seropositive immunity: (days)
P.res = 0.2;                  % (20) seronegative proportion (resistance)

% clinical parameters
%--------------------------------------------------------------------------
P.Tic = 4;                    % (21) asymptomatic period (days)
P.Tsy = 5;                    % (22) symptomatic period  (days)
P.Trd = 16;                   % (23) severe symptoms     (days)

P.sev = 0.01;                 % (24) P(ARDS | symptoms)
P.lat = 1;                    % (25) changes in severity
P.fat = 0.5;                  % (26) P(fatality | ARDS)
P.sur = 1;                    % (27) changes in severity

% testing parameters
%--------------------------------------------------------------------------
P.ttt = 0.036;                % (28) FTTI efficacy
P.tes = [16 8];               % (29) bias (for infection): PCR (Pill. 1 & 2)
P.tts = 1;                    % (30) bias (for infection): LFD
P.del = 3;                    % (31) test delay (days)
P.vac = 32;                   % (32) seroconversion: vaccine (days)
P.fnr = [0.08 0.06];          % (33) false-negative rate  (infected/ious]
P.fpr = [0.0002 0.003];       % (34) false-positive rate: (Sus. and Ab +ve)

P.lim = [1/2 1   2   2]/1000; % (35) testing: capacity (P1, P2, & LFD)
P.rat = [8   24  8   4];      % (36) testing: dispersion
P.ons = [100 200 300 400];    % (37) testing: onset

P.lag = [1 1];                % (38) reporting lag
P.inn = 1;                    % (39) seasonal phase
P.mem = [128 32];             % (40) vaccination (days)
P.rol = [exp(-16) 32];        % (41) vaccination rollout (1st)
P.fol = [exp(-16) 32];        % (42) vaccination rollout (2nd)

P.vef = 0.4;                  % (43) 1 - vaccine efficacy: infection
P.lnk = 0.2;                  % (44) 1 - vaccine efficacy: pathogenicity
P.ves = 0.1;                  % (45) 1 - vaccine efficacy: transmission
P.lnf = 0.05;                 % (46) 1 - vaccine efficacy: fatality

P.con = 0.2;                  % (47) LFD confirmation
P.iso = 8.0;                  % (48) self-isolation (days)
P.Tnn = 256;                  % (49) loss of T-cell immunity (days)

i     = ceil((datenum(date) - datenum('01-02-2020'))/64);

P.lnr = 0.46;                 % (50) LFD sensitivity
P.lpr = 0.0002;               % (51) LFD specificity
P.rel = 1;                    % (52) PCR testing of fatalities
P.pro = 1;                    % (53) contact rate decay (days)
P.oth = 0.1;                  % (54) relative survival outside hospital
P.iad = 1;                    % (55) exponent: transfer to CCU
P.tra = ones(1,i)*64/800;     % (56) transmissibility parameters
P.dps = [2 1];                % (57) doses per seroconversion
P.abs = 1;                    % (58) age-related testing
P.iss = 1;                    % (59) probability of self-isolation
P.rut = 1;                    % (60) Sensitivity:contact rate fluctuations


% infection fatality (for susceptible population)
%--------------------------------------------------------------------------
% IFR (hospital): P.sev*P.fat*100
% IFR (hospital): exp(Ep.sev)*exp(Ep.fat)*100

% Variances (mildly informative priors, apart from initial cases and size)
%==========================================================================
U     = exp( 2);              % flat priors
V     = exp(-2);              % uninformative priors
W     = exp(-4);              % informative priors (weak)
X     = exp(-6);              % informative priors (strong)
Z     = exp(-8);              % informative priors (very strong)

C.N   = U;                    % (01) population size (millions)
C.n   = U;                    % (02) initial cases (cases)
C.r   = X;                    % (03) pre-existing immunity (proportion)
C.o   = W;                    % (04) initial exposed proportion
C.m   = W;                    % (05) relative eflux

% location (exposure) parameters
%--------------------------------------------------------------------------
C.out = X;                    % (06) P(leaving home)
C.sde = W;                    % (07) time constant of lockdown
C.qua = X;                    % (08) time constant of unlocking
C.exp = W;                    % (09) viral spreading (rate)
C.hos = W;                    % (10) admission rate (hospital)
C.ccu = W;                    % (11) admission rate (CCU)
C.s   = X;                    % (12) changes in infectivity

% infection (transmission) parameters
%--------------------------------------------------------------------------
C.Nin = X;                    % (13) effective number of contacts: home
C.Nou = X;                    % (14) effective number of contacts: work
C.trn = X;                    % (16) transmission strength
C.trm = W;                    % (15) seasonality
C.Tin = Z;                    % (17) infected period (days)
C.Tcn = Z;                    % (18) infectious period (days)
C.Tim = X;                    % (19) seropositive immunity (days)
C.res = X;                    % (20) seronegative immunity (proportion)

% clinical parameters
%--------------------------------------------------------------------------
C.Tic = X;                    % (21) asymptomatic period (days)
C.Tsy = X;                    % (22) symptomatic period  (days)
C.Trd = W;                    % (23) CCU period (days)
C.sev = V;                    % (24) P(ARDS | symptoms)
C.lat = V;                    % (25) changes in severity
C.fat = V;                    % (26) P(fatality | ARDS)
C.sur = V;                    % (27) changes in severity

% testing parameters
%--------------------------------------------------------------------------
C.ttt = W;                    % (28) FTTI efficacy
C.tes = V;                    % (29) testing: bias (early)
C.tts = V;                    % (30) testing: bias (late)
C.del = X;                    % (31) test delay (days)
C.vac = X;                    % (32) vaccination time constant (days)
C.fnr = X;                    % (33) false-negative rate
C.fpr = X;                    % (34) false-positive rate

C.lim = V;                    % (35) testing: capacity
C.rat = X;                    % (36) testing: constant (days)
C.ons = U;                    % (37) testing: onset (days)

C.lag = V;                    % (38) reporting lag
C.inn = V;                    % (39) seasonal phase
C.mem = Z;                    % (40) vaccination (initial doses)
C.rol = X;                    % (41) vaccination rollout (1st)
C.fol = X;                    % (42) vaccination rollout (2nd)

C.vef = X;                    % (43) vaccine efficacy: infection
C.lnk = X;                    % (44) vaccine efficacy: pathogenicity
C.ves = X;                    % (45) vaccine efficacy: transmission
C.lnf = X;                    % (46) vaccine efficacy: fatality

C.con = V;                    % (47) LFD confirmation
C.iso = Z;                    % (48) self-isolation (days)
C.Tnn = Z;                    % (49) loss of T-cell immunity (days)

C.lnr = X;                    % (50) LFD sensitivity
C.lpr = X;                    % (51) LFD specificity
C.rel = W;                    % (52) PCR testing of fatalities
C.pro = W;                    % (53) contact rate decay (days)
C.oth = W;                    % (54) relative survival outside hospital
C.iad = V;                    % (55) exponent: transfer to CCU
C.tra = X;                    % (56) transmissibility parameters
C.dps = X;                    % (57) doses per seroconversion
C.abs = X;                    % (58) age-related testing
C.iss = X;                    % (59) probability of self-isolation
C.rut = X;                    % (60) Sensitivity:contact rate fluctuations


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