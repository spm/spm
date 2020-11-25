function [DCM] = DEM_COVID_LTLA
% FORMAT [DCM] = DEM_COVID_LTLA
%
% Demonstration of COVID-19 modelling
%__________________________________________________________________________
%
% This demonstration routine fixed multiple regional death by date and new
% cases data and compiles estimates of latent states for local
% authorities.
%
% Technical details about the dynamic causal model used here can be found
% at https://www.fil.ion.ucl.ac.uk/spm/covid-19/.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DEM_COVID_UTLA.m 8005 2020-11-06 19:37:18Z karl $



% load (ONS) testing death-by-date data
%==========================================================================
url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&metric=newDeaths28DaysByDeathDate&format=csv';
U   = webread(url);
P   = importdata('LADCodesPopulation2019.xlsx');

% get population by lower tier local authority
%--------------------------------------------------------------------------
PN       = P.data(2:end,1);
PCode    = P.textdata(2:end,1);
PName    = P.textdata(:,2);

% get new cases by (lower tier) local authority
%--------------------------------------------------------------------------
AreaCode = table2cell(U(:,3));
AreaType = table2cell(U(:,2));
AreaDate = table2cell(U(:,1));
AreaName = table2cell(U(:,4));

j        = find(ismember(AreaType,'ltla') & ~ismember(AreaCode,'null') );
AreaCode = AreaCode(j);
AreaDate = AreaDate(j);
AreaName = AreaName(j);
AreaCase = table2array(U(j,5));
AreaMort = table2array(U(j,6));

% organise via NHS trust
%--------------------------------------------------------------------------
Area  = unique(AreaCode);
k     = 0;
for i = 1:numel(Area)
    
    % get code
    %----------------------------------------------------------------------
    j = find(ismember(AreaCode,Area(i)));
    
    if any(isfinite(AreaCase(j))) && any(isfinite(AreaMort(j)))
        k           = k + 1;
        D(k).cases  = AreaCase(j);
        D(k).deaths = AreaMort(j);
        D(k).date   = datenum([AreaDate{j}]);
        D(k).name   = AreaName(j(1));
        D(k).code   = Area(i);
        D(k).N      = PN(find(ismember(PCode,Area(i)),1));
        
        if isempty(D(k).N)
            k = k - 1; disp(Area(i));
        end
    end
    
end

% clear
%--------------------------------------------------------------------------
clear P PN PCD AreaCase AreaCode AreaDate AreaMort AreaName AreaType

% dates to generate
%--------------------------------------------------------------------------
d0     = min(spm_vec(D.date));
d0     = min(d0,datenum('01-Feb-2020','dd-mmm-yyyy'));
M.date = datestr(d0,'dd-mmm-yyyy');
dates  = d0:max(spm_vec(D.date));
    
    
% fit each regional dataset
%==========================================================================
for r = 1:numel(D)
    
    fprintf('%d out of %d\n',r,numel(D));
    try, clear Y; end
    
    % created data structure
    %----------------------------------------------------------------------
    Y(1).type = 'PCR cases (ONS)'; % daily PCR positive cases (by specimen)
    Y(1).unit = 'number/day';
    Y(1).U    = 2;
    Y(1).date = D(r).date;
    Y(1).Y    = D(r).cases;
    Y(1).h    = 0;
    
    Y(2).type = 'Daily deaths (ONS: 28-days)'; % daily covid-related deaths (28 days)
    Y(2).unit = 'number/day';
    Y(2).U    = 1;
    Y(2).date = D(r).date;
    Y(2).Y    = D(r).deaths;
    Y(2).h    = 2;

    % remove NANs and sort by date
    %----------------------------------------------------------------------
    [Y,YS] = spm_COVID_Y(Y,M.date);
    
    % data structure with vectorised data and covariance components
    %----------------------------------------------------------------------
    xY.y  = spm_vec(Y.Y);
    xY.Q  = spm_Ce([Y.n]);
    hE    = spm_vec(Y.h);
    
    % get and set priors
    %----------------------------------------------------------------------
    [pE,pC] = spm_SARS_priors;
   
    pE.N   = log(D(r).N/1e6);      % population of local authority
    pC.N   = 0;
    pE.n   = -8;                   % initial number of cases (n)

    % model specification
    %======================================================================
    M.Nmax = 16;                   % maximum number of iterations
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
    DCM(r).F  = F;
    
    % now-casting for this region and date
    %======================================================================
    H     = spm_figure('GetWin',D(r).name{1}); clf;
    %----------------------------------------------------------------------
    M.T   = numel(dates) + 32;
    [Y,X] = spm_SARS_gen(DCM(r).Ep,M,[1 2]);
    spm_SARS_plot(Y,X,YS);
    
    
    %----------------------------------------------------------------------
    % Y(:,1)  - number of new deaths
    % Y(:,2)  - number of new cases
    % Y(:,3)  - CCU bed occupancy
    % Y(:,4)  - effective reproduction rate (R)
    % Y(:,5)  - population immunity (%)
    % Y(:,6)  - total number of tests
    % Y(:,7)  - contagion risk (%)
    % Y(:,8)  - prevalence of infection (%)
    % Y(:,9)  - number of infected at home, untested and asymptomatic
    % Y(:,10) - new cases per day
    %----------------------------------------------------------------------
    M.T     = numel(dates);
    Y       = spm_SARS_gen(DCM(r).Ep,M,[4 5 8 9 10]);
    
    DR(:,r) = Y(:,1);                       % Reproduction ratio
    DI(:,r) = Y(:,2);                       % Prevalence of immunity
    DP(:,r) = Y(:,3);                       % Prevalence of infection
    DC(:,r) = Y(:,4);                       % Infected, asymptomatic people
    DT(:,r) = Y(:,5)/exp(Ep.N)/10;          % New daily cases per 100K
    
    
    % supplement with table of posterior expectations
    %----------------------------------------------------------------------
    subplot(3,2,2), cla reset, axis([0 1 0 1])
    title(D(r).name,'Fontsize',16)
    
    str = sprintf('Population: %.2f million',exp(Ep.N));
    text(0,0.9,str,'FontSize',10,'FontWeight','bold','Color','k')
    
    str = sprintf('Reproduction ratio: %.2f',DR(end,r));
    text(0,0.8,str,'FontSize',10,'FontWeight','bold','Color','k')
    
    str = sprintf('Infected, asymptomatic people: %.0f',DC(end,r));
    text(0,0.7,str,'FontSize',10,'FontWeight','bold','Color','k')
    
    str = sprintf('Daily new cases: %.0f per 100,000',DT(end,r));
    text(0,0.6,str,'FontSize',10,'FontWeight','bold','Color','r')
    
    str = sprintf('Prevalence of infection: %.2f%s',DP(end,r),'%');
    text(0,0.5,str,'FontSize',10,'FontWeight','bold','Color','r')
    
    str = sprintf('Prevalence of immunity: %.1f%s',DI(end,r),'%');
    text(0,0.4,str,'FontSize',10,'FontWeight','bold','Color','k')
    
    str = {'The prevalences refer to the estimated population ' ...
           '(based on ONS census figures for lower tier local authorities)'};
    text(0,0.0,str,'FontSize',8,'Color','k')
    
    spm_axis tight, axis off
    drawnow
    savefig(H,[strrep(strrep(strrep(D(r).name{1},' ','_'),',',''),'''',''),'.fig']);
    close(H);
    
end

% save
%----------------------------------------------------------------------
try, clear ans, end
try, clear H,   end
save COVID_LA



