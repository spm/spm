function [DCM] = DEM_COVID_AGE
% FORMAT [DCM] = DEM_COVID_AGE
% LA - local authority
%
% Demonstration of COVID-19 modelling
%__________________________________________________________________________
%
% This demonstration routine fits multiple regional death by date and new
% cases data and compiles estimates of latent states for local
% authorities.
%
% Technical details about the dynamic causal model used here can be found
% at https://www.fil.ion.ucl.ac.uk/spm/covid-19/.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

%% age demographics
%--------------------------------------------------------------------------
sage = {'00_04','05_09','10_14','15_19','20_24','25_29','30_34','35_39',...
        '40_44','45_49','50_54','55_59','60_64','65_69','70_74','75_79', ...
        '80_84','85_89','90+'};

page = [3.86, 4.15, 3.95, 3.66, 4.15, 4.51, 4.50, 4.40, 4.02, 4.4, 4.66, ...
        4.41, 3.76, 3.37, 3.32, 2.33, 1.72, 1.04, 0.61];

EnglandUK  = 66.79/56.28;


% download from web options
%--------------------------------------------------------------------------
options   = weboptions('ContentType','table');
options.Timeout = 20;

% load (ONS) testing death-by-date data
%==========================================================================
url       = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv';
deaths    = webread(url,options);
url       = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newCasesBySpecimenDateAgeDemographics&format=csv';
cases     = webread(url,options);

% get new cases by (lower tier) local authority
%--------------------------------------------------------------------------
DeathDate = datenum(table2array(deaths(:,1)));
DeathAge  = table2array(deaths(:,6));
DeathNum  = table2array(deaths(:,7));

CaseDate  = datenum(table2array(cases(:,1)));
CaseAge   = table2array(cases(:,6));
CaseNum   = table2array(cases(:,7));


% get empirical priors from national data
%==========================================================================
if false
    DCM = DEM_COVID_UK;
    save('DCM_UK.mat','DCM')
    PCM = DCM;
else
    PCM = load('DCM_UK.mat','DCM');
    PCM = PCM.DCM;
end

% dates to generate
%--------------------------------------------------------------------------
d0     = min(spm_vec([DeathDate; CaseDate]));
d0     = min(d0,datenum(PCM.M.date,'dd-mm-yyyy'));
M.date = datestr(d0,'dd-mm-yyyy');
dates  = d0:max(spm_vec([DeathDate; CaseDate]));

% free parameters of local model (fixed effects)
%==========================================================================
free  = {'n','res','sev','lat','fat','tes','tts','pcr'};

% (empirical) prior expectation
%--------------------------------------------------------------------------
pE    = PCM.Ep;

% (empirical) prior covariances
%--------------------------------------------------------------------------
pC    = spm_zeros(PCM.M.pC);
for i = 1:numel(free)
    pC.(free{i}) = spm_zeros(PCM.M.pC.(free{i})) + 1;
end

% fit each age bound
%==========================================================================
ages    = unique(sage);
age{1}  = ages(1:4);
age{2}  = ages(5:9);
age{3}  = ages(10:13);
age{4}  = ages(14:end);

Nage(1) = sum(page(1:4));
Nage(2) = sum(page(5:9));
Nage(3) = sum(page(10:13));
Nage(4) = sum(page(14:end));

for r = 1:numel(age)
    
    fprintf('%d out of %d\n',r,numel(age));
    try, clear Y; end
    
    % get age bands
    %----------------------------------------------------------------------
    D     = 0;
    for i = 1:numel(age{r})
        j = find(ismember(CaseAge,age{r}{i}));
        D = D + CaseNum(j);
    end
    
    % created data structure
    %----------------------------------------------------------------------
    Y(1).type = 'PCR cases (ONS)'; % daily PCR positive cases (by specimen)
    Y(1).unit = 'number/day';
    Y(1).U    = 2;
    Y(1).date = CaseDate(j);
    Y(1).Y    = D;
    Y(1).h    = 0;
    
    % get age bands
    %----------------------------------------------------------------------
    D     = 0;
    for i = 1:numel(age{r})
        j = find(ismember(DeathAge,age{r}{i}));
        D = D + DeathNum(j);
    end
    
    Y(2).type = 'Daily deaths (ONS: 28-days)'; % daily covid-related deaths (28 days)
    Y(2).unit = 'number/day';
    Y(2).U    = 15;
    Y(2).date = DeathDate(j);
    Y(2).Y    = D*EnglandUK;
    Y(2).h    = 4;
    
    % remove NANs, smooth and sort by date
    %----------------------------------------------------------------------
    [Y,S]  = spm_COVID_Y(Y,M.date,16);
    
    % data structure with vectorised data and covariance components
    %----------------------------------------------------------------------
    xY.y   = spm_vec(Y.Y);
    xY.Q   = spm_Ce([Y.n]);
    hE     = spm_vec(Y.h);
    
    % get and set priors
    %----------------------------------------------------------------------
    pE.N   = log(Nage(r));         % population of age band
    
    % model specification
    %======================================================================
    M.Nmax = 32;                   % maximum number of iterations
    M.G    = @spm_SARS_gen;        % generative function
    M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
    M.pE   = pE;                   % prior expectations (parameters)
    M.pC   = pC;                   % prior covariances  (parameters)
    M.hE   = hE;                   % prior expectation  (log-precision)
    M.hC   = 1/512;                % prior covariances  (log-precision)
    M.T    = Y;                    % data structure
    U      = spm_vec(Y.U);         % outputs to model
    
    % model inversion with Variational Laplace (Gauss Newton)
    %----------------------------------------------------------------------
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,xY);
    
    % save prior and posterior estimates (and log evidence)
    %----------------------------------------------------------------------
    DCM(r).M  = M;
    DCM(r).Ep = Ep;
    DCM(r).Eh = Eh;
    DCM(r).Cp = Cp;
    DCM(r).Y  = Y;
    DCM(r).xY = xY;
    DCM(r).S  = S;
    DCM(r).F  = F;
    
    % now-casting for this region and date
    %======================================================================
    H      = spm_figure('GetWin',['Age range: ' age{r}{1}]); clf;
    %----------------------------------------------------------------------
    M.T    = numel(dates) + 32;
    u      = [find(U == 1)];
    S      = DCM(r).S(:,u); S(isnan(S)) = 0;
    [Y,X]  = spm_SARS_gen(DCM(r).Ep,M,1);
    spm_SARS_plot(Y,X,S,1);
    
end

%% parametric intervention (vaccination)
%==========================================================================
spm_figure('GetWin','testing and cases'); clf; subplot(2,1,1)
period = {M.date,'01-06-2021'};         % duration of epidemic

% targeted vaccination rate
%--------------------------------------------------------------------------
rol           = exp(Ep.rol(1))*sum(Nage)/Nage(end);


% no vaccine rollout
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = 'rol';
NPI(1).Q      = [1e-8 exp(Ep.rol(2:3))];
NPI(1).dates  = {'20-12-2020',period{2}};

% vaccine rollout based upon posterior estimates
%--------------------------------------------------------------------------
NPI(2).period = period;
NPI(2).param  = 'rol';
NPI(2).Q      = [rol exp(Ep.rol(2:3))];
NPI(2).dates  = {'20-12-2020',period{2}};


% solve for both scenarios
%--------------------------------------------------------------------------
M.T   = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
T     = (1:M.T) + dates(1) - 1;
for i = 1:2
    
    YY    = 0;                               % predicted cumulative deaths
    SS    = 0;                               % empirical cumulative deaths
    VV    = 0;                               % cumulative vaccinations
    for r = 1:numel(DCM)
        
        % limit vaccination to oldest age group
        %------------------------------------------------------------------
        if r == numel(DCM)
            npi = NPI(i);
        else
            npi = NPI(1);
        end
        
        % unpack model and posterior expectations
        %------------------------------------------------------------------
        Ep  = DCM(r).Ep;                            % posterior expectation
        S   = DCM(r).S;                             % smooth timeseries
        
        % quantify the effect of efficient intervention (NPI)
        %==================================================================
        [Y,X] = spm_SARS_gen(Ep,M,[15 22],npi);
        
        % record number of deaths, vaccinations and empirical mortality
        %------------------------------------------------------------------
        YY  = YY + Y(:,1);
        VV  = VV + Y(:,2);
        SS  = SS + S(:,2);
        
        % plot for vaccination scenario
        %------------------------------------------------------------------
        if i > 1
            fill([T(:); flipud(T(:))],[0*YY(:); flipud(YY(:))],'r','EdgeColor','none','FaceAlpha',0.4)
            if r == numel(DCM)
                plot(dates,SS,'.','Color','k'), hold on
            end
        end
        
        % add counterfactual (no vaccination)
        %------------------------------------------------------------------
        if r == numel(DCM) && (i == 1)
            plot(T,YY,'LineStyle','-.','Color','r'), hold on
        end
        datetick('x','dd-mmm','keeplimits','keepticks')
        
    end
    
    % add text and legend
    %------------------------------------------------------------------
    if i > 1
        d    = datenum('01-Mar-2021','dd-mmm-yyyy');
        CV   = VV;
        DV   = diff(CV)*1e6;
        DV   = DV(DV > 1);
        j    = find(T >= d,1);
        
        str  = sprintf('%2.1f M by %s',CV(j),datestr(T(j),'dd-mmm'));
        leg  = sprintf('%2.0f per day (max: %2.0f)',mean(DV),max(DV));
        
        CV   = 100*CV;
        plot(T,CV,'b'), hold on
        text(T(j),CV(j),str)
        title(leg), xlabel('date'), ylabel('daily deaths (weekly average)')
        box off
    end
    
end

% ancillary figures
%--------------------------------------------------------------------------
d    = datenum('31-Jan-2021','dd-mmm-yyyy');
disp(YY(find(T >= d,1)))
d    = datenum('15-Feb-2021','dd-mmm-yyyy');
disp(YY(find(T >= d,1)))
d    = datenum('01-Mar-2021','dd-mmm-yyyy');
disp(YY(find(T >= d,1)))

sprintf('%2.0f per week (max: %2.0f)',mean(DV)*7,max(DV)*7)

% Notes for Channel 4
%==========================================================================
vaccination_rate = 2431648 - 1296432;       % vaccination rate per week
datestr(datenum('20-Dec-2020','dd-mmm-yyyy') + round(13.9e6/(vaccination_rate/7)))





