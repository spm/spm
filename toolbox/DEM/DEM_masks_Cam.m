function DEM_masks_Cam
% Bayesian model reduction for COVID models
% FORMAT DEM_mask_Cam
%
% This subroutine applies Bayesian model reduction to a DCM for COVID,
% asking whether mask wearing can be treated as fixed parameters by
% reducing its prior variance to 0.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% setup
%==========================================================================

% get posterior estimates
%--------------------------------------------------------------------------
clear
cd('C:\Users\karl\Dropbox\Coronavirus\Dashboard')
DCM = load('DCM_UK_tmp.mat','DCM');
DCM = DCM.DCM;

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
pE  = DCM.M.pE;                              % prior expectation
pC  = DCM.M.pC;                              % prior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs
A   = DCM.A;                                 % age cohort

Vp  = spm_unvec(spm_vec(diag(Cp)),Ep);       % posterior variances

% Bayesian model comparison
%==========================================================================
param = {'msk','tra'};

qE    = Ep;
qC    = Cp;
L     = [];
str   = {};
nP    = numel(param);
for i = 1:nP

    % can we ignore this effect? (i.e., vanishing prior variance)
    %------------------------------------------------------------------
    rE    = pE;
    rC    = pC;
    for j = 1:numel(pE.(param{i}))
        rE.(param{i})(j) = rE.(param{i})(j) - 2;
        rC.(param{i})(j) = rC.(param{i})(j)/exp(8);
        F  = spm_log_evidence(qE,qC,pE,pC,rE,rC);
        L(end + 1)   = 0 - F;
        str{end + 1} = param{i};
    end
end

% plot and display results of Bayesian model comparison
%--------------------------------------------------------------------------
spm_figure('GetWin','BMR-MASKS'); clf;
subplot(2,1,1)
bar(L), hold on
plot([0 numel(L)],[5,5],'-.')
set(gca,'XTickLabel',str,'YLim',[-1024 1024])
title('Bayesian model comparison','Fontsize',16)
xlabel('parameter'),ylabel('log Bayes factor')
axis square, box off

subplot(2,2,3), hold off
x     = linspace(exp(-8),8/100,64);
for i = 1:numel(Ep.msk)
    p = spm_Npdf(log(x),Ep.msk(i),Vp.msk(i));
    plot(x*100,p), hold on
end
title('Efficacy of face covering','Fontsize',16)
xlabel('efficacy: % reduction in tranmission'), ylabel('posterior probability')
xlabel('DCM parameter'), ylabel('correlation')
axis square, box off, legend({'children','young','adults','elderly'})

subplot(2,2,4), hold off
j     = spm_fieldindices(pE,'msk');
Rp    = spm_cov2corr(Cp);
for i = j(:)'
    plot(Rp(:,i)), hold on
end
title('Posterior correlations','Fontsize',16)
xlabel('DCM parameter'), ylabel('correlation')
axis square, box off


spm_fieldindices(pE,65)
spm_fieldindices(pE,148)

% prior expectation
%==========================================================================
% P.N   = 64;                   % (01) population size (millions)
% P.n   = exp(-8);              % (02) proportion of initial cases (cases)
% P.r   = 0.1;                  % (03) pre-existing immunity (proportion)
% P.o   = 0.1;                  % (04) initial exposed proportion
% P.m   = 0.1;                  % (05) relative eflux
% 
% % location (exposure) parameters
% %--------------------------------------------------------------------------
% P.out = 0.4;                  % (06) P(leaving home)
% P.sde = 16;                   % (07) time constant of lockdown
% P.qua = 128;                  % (08) time constant of unlocking
% P.exp = 0.02;                 % (09) viral spreading (rate)
% P.hos = 1;                    % (10) admission rate (hospital) [erf]
% P.ccu = 0.1;                  % (11) admission rate (CCU)
% P.s   = [1 1 1 1 1 1];        % (12) infectivity exponents
% 
% % infection (transmission) parameters
% %--------------------------------------------------------------------------
% P.Nin = 2;                    % (13) effective number of contacts: home
% P.Nou = 24;                   % (14) effective number of contacts: work
% P.trn = 0.2;                  % (15) transmission strength (secondary attack rate)
% P.trm = 0.04;                 % (16) seasonality
% P.Tin = 3;                    % (17) infected period (days)
% P.Tcn = 4;                    % (18) infectious period (days)
% P.Tim = 128;                  % (19) seropositive immunity: (days)
% P.res = 0.2;                  % (20) seronegative proportion (resistance)
% 
% % clinical parameters
% %--------------------------------------------------------------------------
% P.Tic = 4;                    % (21) asymptomatic period (days)
% P.Tsy = 5;                    % (22) symptomatic period  (days)
% P.Trd = 16;                   % (23) severe symptoms     (days)
% 
% P.sev = 0.01;                 % (24) P(ARDS | symptoms)
% P.lat = 1;                    % (25) changes in severity
% P.fat = 0.5;                  % (26) P(fatality | ARDS)
% P.sur = 1;                    % (27) changes in severity
% 
% % testing parameters
% %--------------------------------------------------------------------------
% P.ttt = 0.036;                % (28) FTTI efficacy
% P.tes = [16 8];               % (29) bias (for infection): PCR (Pill. 1 & 2)
% P.tts = 1;                    % (30) bias (for infection): LFD
% P.del = 3;                    % (31) test delay (days)
% P.vac = 32;                   % (32) seroconversion: vaccine (days)
% P.fnr = [0.08 0.06];          % (33) false-negative rate  (infected/ious]
% P.fpr = [0.0002 0.003];       % (34) false-positive rate: (Sus. and Ab +ve)
% 
% P.lim = [1/2 1   2   2]/1000; % (35) testing: capacity (P1, P2, & LFD)
% P.rat = [8   24  8   4];      % (36) testing: dispersion
% P.ons = [100 200 300 400];    % (37) testing: onset
% 
% P.lag = [1 1];                % (38) reporting lag
% P.inn = 1;                    % (39) seasonal phase
% P.mem = [128 32];             % (40) vaccination (days)
% P.rol = [exp(-16) 32];        % (41) vaccination rollout (1st)
% P.fol = [exp(-16) 32];        % (42) vaccination rollout (2nd)
% 
% P.vef = 0.4;                  % (43) 1 - vaccine efficacy: infection
% P.lnk = 0.2;                  % (44) 1 - vaccine efficacy: pathogenicity
% P.ves = 0.1;                  % (45) 1 - vaccine efficacy: transmission
% P.lnf = 0.05;                 % (46) 1 - vaccine efficacy: fatality
% 
% P.con = 0.2;                  % (47) LFD confirmation
% P.iso = 8.0;                  % (48) self-isolation (days)
% P.Tnn = 256;                  % (49) loss of T-cell immunity (days)
% P.lnr = 0.46;                 % (50) LFD sensitivity
% P.lpr = 0.0002;               % (51) LFD specificity
% P.rel = 1;                    % (52) PCR testing of fatalities
% P.pro = 1;                    % (53) contact rate decay (days)
% P.oth = 0.1;                  % (54) relative survival outside hospital
% P.iad = 1;                    % (55) exponent: transfer to CCU
% P.tra = ones(1,nb)/8;         % (56) transmissibility parameters
% P.dps = [2 1];                % (57) doses per seroconversion
% P.abs = 1;                    % (58) age-related testing
% P.iss = 1;                    % (59) probability of self-isolation
% P.rut = 1;                    % (60) Sensitivity:contact rate fluctuations
% P.msk = 0.01;                 % (61) Sensitivity:contact rate fluctuations
%--------------------------------------------------------------------------

%% Outcomes
%==========================================================================
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
% Y(:,29) - Proportion vaccinated
% Y(:,30) - Cumulative admissions
% Y(:,31) - Vaccine effectiveness (prevalence)
% Y(:,32) - Gross domestic product

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
global CHOLD GHOLD
CHOLD  = 1;
GHOLD  = 1;
%--------------------------------------------------------------------------
M.T    = datenum(M.end) - datenum(M.date);
u      = 1;                                  % empirical outcome
a      = 0;                                  % age cohort (0 for everyone)

Ep.msk = DCM.Ep.msk + 0;                     % adjusted (log) parameter
[Z ,X] = spm_SARS_gen(Ep,M,u,[],a);          % posterior prediction

% and plot
%--------------------------------------------------------------------------
j     = [];
for i = 1:numel(u)
    j = [j find(U == u(i) & A == a(i))];
end
spm_SARS_plot(Z,X,S(:,j),u)

% reduce mask wearing
%--------------------------------------------------------------------------
Ep.msk = DCM.Ep.msk - 4;                     % adjusted (log) parameter
[Z0,X] = spm_SARS_gen(Ep,M,u,[],a);          % posterior prediction

% and plot
%--------------------------------------------------------------------------
spm_SARS_plot(Z0,X,S(:,j),u)

% overlay mask wearing
%--------------------------------------------------------------------------
mask  = [];
for i = datenum(M.date):7:datenum(M.end)
    mask(end + 1) = DEM_COVID_MASKS(i);
end
subplot(3,2,1)
plot(max(Z0)*mask/max(mask))

disp('cases saved:')
disp(sum(Z0) - sum(Z))
disp('percent:')
disp(100*(sum(Z0) - sum(Z))/sum(Z))

% repeat with confidence (credible) intervals
%==========================================================================
spm_figure('GetWin','CI'); clf;
subplot(2,1,1)
Ep.msk = DCM.Ep.msk + 0;                     % adjusted (log) parameter
spm_SARS_ci(Ep,Cp,S(1:M.T,j),u,M); hold on
Ep.msk = DCM.Ep.msk - 4;                     % adjusted (log) parameter
spm_SARS_ci(Ep,Cp,[],u,M); hold on

% overlay mask wearing
%--------------------------------------------------------------------------
mask  = [];
t     = datenum(M.date):datenum(M.end);
for i = t
    mask(end + 1) = DEM_COVID_MASKS(i);
end
plot(t,max(Z0)*mask/max(mask))








