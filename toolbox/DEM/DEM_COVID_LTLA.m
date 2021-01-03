function [DCM] = DEM_COVID_LTLA(LA)
% FORMAT [DCM] = DEM_COVID_LTLA(LA)
% LA - local authority
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

% download from web options
%--------------------------------------------------------------------------
options = weboptions('ContentType','table');
options.Timeout = 20;

% load (ONS) testing death-by-date data
%==========================================================================
url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newDeaths28DaysByDeathDateRate&metric=newOnsDeathsByRegistrationDate&metric=uniqueCasePositivityBySpecimenDateRollingSum&metric=uniquePeopleTestedBySpecimenDateRollingSum&metric=newCasesBySpecimenDate&format=csv';
U   = webread(url,options);
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
newDeath = table2array(U(j,5));
OnsDeath = table2array(U(j,6));
Positive = table2array(U(j,7));
Tested   = table2array(U(j,8));
newCases = table2array(U(j,9));

% organise via NHS trust
%--------------------------------------------------------------------------
Area  = unique(AreaCode);
k     = 0;
for i = 1:numel(Area)
    
    % get code
    %----------------------------------------------------------------------
    j = find(ismember(AreaCode,Area(i)));
    
    if any(isfinite(newCases(j))) && (any(isfinite(newDeath(j))) || any(isfinite(OnsDeath(j))))
        k           = k + 1;
        D(k).cases  = newCases(j);
        D(k).deaths = newDeath(j);
        D(k).cert   = OnsDeath(j);
        D(k).tests  = Tested(j);
        D(k).rate   = Positive(j);
        
        D(k).date   = datenum([AreaDate{j}]);
        D(k).name   = AreaName(j(1));
        D(k).code   = Area(i);
        D(k).N      = PN(find(ismember(PCode,Area(i)),1));
        
        % replace early NaNs with zero
        %------------------------------------------------------------------
        j = D(k).date(:) < mean(D(k).date) & isnan(D(k).cases);
        D(k).cases(j)  = 0;
        j = D(k).date(:) < mean(D(k).date) & isnan(D(k).deaths);
        D(k).deaths(j) = 0;
        
        if isempty(D(k).N)
            k = k - 1; disp(Area(i));
        end
    end
    
end

% clear
%--------------------------------------------------------------------------
clear P PN PCD AreaCase AreaCode AreaDate AreaMort AreaName AreaType

% check requested if necessary
%--------------------------------------------------------------------------
if nargin
    try
        i = find(ismember([D.name],LA));
    catch
        i = find(ismember([D.code],LA));
    end
    if isempty(i)
        disp('local authority not found, please check code or name')
        return
    else
        D = D(i);
    end
end

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
d0     = min(spm_vec(D.date));
d0     = min(d0,datenum(PCM.M.date,'dd-mmm-yyyy'));
M.date = datestr(d0,'dd-mmm-yyyy');
dates  = d0:max(spm_vec(D.date));

% free parameters of local model (fixed effects)
%==========================================================================
free  = {'n','r','o','m','sde','qua','exp','s','Nin','Nou','tes','tts','lim','ons'};


% (empirical) prior expectation
%--------------------------------------------------------------------------
pE    = PCM.Ep;

% (empirical) prior covariances
%--------------------------------------------------------------------------
pC    = spm_zeros(PCM.M.pC);
for i = 1:numel(free)
    pC.(free{i}) = PCM.M.pC.(free{i});
end

%%% try, D   = D(1:8); end %%%%

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
    Y(1).h    = 2;
    
    % deal with Welsh death data
    %----------------------------------------------------------------------
    if ismember('W',D(r).code{1})
        
        Y(2).type = 'Certified deaths (ONS)'; % weekly covid related deaths
        Y(2).unit = 'number/day';
        Y(2).U    = 1;
        Y(2).date = D(r).date;
        Y(2).Y    = D(r).cert/7;
        Y(2).h    = 0;
    else
        Y(2).type = 'Daily deaths (ONS: 28-days)'; % daily covid-related deaths (28 days)
        Y(2).unit = 'number/day';
        Y(2).U    = 1;
        Y(2).date = D(r).date;
        Y(2).Y    = D(r).deaths;
        Y(2).h    = 0;
    end
    
    Y(3).type = 'PCR positivity (GOV)'; % positivity (England)
    Y(3).unit = 'percent';
    Y(3).U    = 23;
    Y(3).date = D(r).date;
    Y(3).Y    = D(r).rate + 1;
    Y(3).h    = 0;
    
    Y(4).type = 'PCR tests (ONS)'; % daily PCR tests performed
    Y(4).unit = 'number/day';
    Y(4).U    = 6;
    Y(4).date = D(r).date;
    Y(4).Y    = D(r).tests;
    Y(4).h    = 0;
    
    % remove NANs, smooth and sort by date
    %----------------------------------------------------------------------
    [Y,S] = spm_COVID_Y(Y([1,2,3]),M.date,16);
        
    % data structure with vectorised data and covariance components
    %----------------------------------------------------------------------
    xY.y  = spm_vec(Y.Y);
    xY.Q  = spm_Ce([Y.n]);
    hE    = spm_vec(Y.h);
    
    % get and set priors
    %----------------------------------------------------------------------
    pE.N   = log(D(r).N/1e6);      % population of local authority
    
    % model specification
    %======================================================================
    M.Nmax = 32;                   % maximum number of iterations
    M.G    = @spm_SARS_gen;        % generative function
    M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
    M.pE   = pE;                   % prior expectations (parameters)
    M.pC   = pC;                   % prior covariances  (parameters)
    M.hE   = hE + 2;               % prior expectation  (log-precision)
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
    H      = spm_figure('GetWin',D(r).name{1}); clf;
    %----------------------------------------------------------------------
    M.T    = numel(dates) + 32;
    u      = [find(U == 1) find(U == 2) find(U == 23)];
    S      =  DCM(r).S(:,u); S(isnan(S)) = 0;
    [Y,X]  = spm_SARS_gen(DCM(r).Ep,M,[1 2 23]);
    spm_SARS_plot(Y,X,S,[1 2 23]);
    
    
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
    M.P     = Ep;
    Y       = spm_SARS_gen(DCM(r).Ep,M,[4 5 8 9 10]);
    DR(:,r) = Y(:,1);                       % Reproduction ratio
    DI(:,r) = Y(:,2);                       % Prevalence of immunity
    DP(:,r) = Y(:,3);                       % Prevalence of infection
    DC(:,r) = Y(:,4);                       % Infected, asymptomatic people
    DT(:,r) = Y(:,5)*1000;                  % New daily cases per 100K
    
    % cumulative incidence of infection (with 90% credible intervals)
    %----------------------------------------------------------------------
    subplot(3,2,2)
    [S,CS]  = spm_SARS_ci(DCM(r).Ep,DCM(r).Cp,[],10,M);
    CS      = 1.6449*sqrt(diag(CS));
    TT(:,r) = S;                            % cumulative incidence (%)
    TU(:,r) = S + CS;                       % upper bound
    TL(:,r) = S - CS;                       % lower bound
    
    % supplement with table of posterior expectations
    %----------------------------------------------------------------------
    subplot(3,2,2), cla reset, axis([0 1 0 1])
    title(D(r).name,'Fontsize',16)
    
    str = sprintf('Population: %.2f million', exp(Ep.N));
    text(0,0.9,str,'FontSize',10,'FontWeight','bold','Color','k')
    
    str = sprintf('Reproduction ratio: %.2f',          DR(end,r));
    text(0,0.8,str,'FontSize',10,'FontWeight','bold','Color','k')
    
    str = sprintf('Infected asymptomatic people: %.0f',DC(end,r));
    text(0,0.7,str,'FontSize',10,'FontWeight','bold','Color','k')
    
    str = sprintf('Daily new cases: %.0f per 100,000', DT(end,r));
    text(0,0.6,str,'FontSize',10,'FontWeight','bold','Color','r')
    
    str = sprintf('Prevalence of infection: %.2f%s',   DP(end,r),'%');
    text(0,0.5,str,'FontSize',10,'FontWeight','bold','Color','r')
    
    str = sprintf('Prevalence of immunity: %.1f%s',    DI(end,r),'%');
    text(0,0.4,str,'FontSize',10,'FontWeight','bold','Color','k')
    
    str = sprintf('Proportion infected: %.1f (%.1f-%.1f)%s',TT(end,r),TL(end,r),TU(end,r),'%');
    text(0,0.3,str,'FontSize',10,'FontWeight','bold','Color','k')
    
    str = {'The prevalences refer to the estimated population ' ...
        '(based on ONS census figures for lower tier local authorities)'};
    text(0,0.0,str,'FontSize',8,'Color','k')
    
    spm_axis tight, axis off
    drawnow
    if numel(D) > 16
        
        % save figure
        %------------------------------------------------------------------
        savefig(H,[strrep(strrep(strrep(D(r).name{1},' ','_'),',',''),'''',''),'.fig']);
        close(H);
        
        % save working variables
        %------------------------------------------------------------------
        try, clear ans,     end
        try, clear H,       end
        try, save COVID_LA, end
        
    end
end


return


%% regional analysis (one)
%==========================================================================

% get regional data and parameter estimates
%--------------------------------------------------------------------------
load COVID_LA
dn    = 64;
M.T   = numel(dates) + dn;
DC    = zeros(M.T,numel(DCM));
DP    = zeros(M.T,numel(DCM));
DM    = zeros(M.T,numel(DCM));
DV    = zeros(M.T,numel(DCM));
DB    = zeros(M.T,numel(DCM));
for r = 1:numel(DCM)
    
    % now-casting for this region and date
    %----------------------------------------------------------------------
    Y      = spm_SARS_gen(DCM(r).Ep,M,[2 8]);
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
    DC(:,r) = Y(:,1);                       % reported cases
    DP(:,r) = Y(:,2);                       % prevalence 
    
    % remove fluctuations in transmission
    %----------------------------------------------------------------------
    Ep      = DCM(r).Ep;
    Ep.vir  = spm_zeros(Ep.vir);
    Y       = spm_SARS_gen(Ep,M,[2 8]);
    DV(:,r) = Y(:,1);                       % prevalence (transmission)
    
    % remove fluctuations in contact rate
    %----------------------------------------------------------------------
    Ep      = DCM(r).Ep;
    Ep.mob  = spm_zeros(Ep.mob);
    Y       = spm_SARS_gen(Ep,M,[2 8]);
    DM(:,r) = Y(:,1);                       % prevalence (contact rate)
    
    % remove fluctuations in both
    %----------------------------------------------------------------------
    Ep      = DCM(r).Ep;
    Ep.vir  = spm_zeros(Ep.vir);
    Ep.mob  = spm_zeros(Ep.mob);
    Y       = spm_SARS_gen(Ep,M,[2 8]);
    DB(:,r) = Y(:,1);                       % prevalence (both)
    
    disp(r)

end
    
%% Rank regions by first and second peaks (before and after 200 days)
%==========================================================================
spm_figure('GetWin','Variation'), clf
%--------------------------------------------------------------------------
n     = numel(dates);
t1    = 1:200;
t2    = 201:n;
YC    = zeros(n,numel(DCM));
YD    = zeros(n,numel(DCM));
for i = 1:numel(DCM)
    
    % vectorised data and parameters
    %----------------------------------------------------------------------
    EP(:,i)  = spm_vec(DCM(i).Ep);
    
    % deal with missing data
    %----------------------------------------------------------------------
    S        = DCM(i).S;
    S(isnan(S)) = 0;
    j        = (1:size(S,1)) + n - size(S,1);
    YC(j,i)  = S(:,1);                     % new cases
    
end

% plot reported cases
%--------------------------------------------------------------------------
YC    = [YC; zeros(dn,size(YC,2))];
n     = size(YC,1);
t1    = 1:200;
t2    = 201:n;
Y1    = zeros(numel(t1),numel(DCM));
Y2    = zeros(numel(t2),numel(DCM));
for i = 1:numel(DCM)
    
    % normalise to peak amplitude
    %----------------------------------------------------------------------
    Y       = YC;
    Y1(:,i) = Y(t1,i)/max(Y(t1,i));
    Y2(:,i) = Y(t2,i)/max(Y(t2,i));
    
    % amplitude and time of peaks
    %----------------------------------------------------------------------
    [m,j]   = max(Y(t1,i));
    T1(i,1) = j;
    
    [m,j]   = max(Y(t2,i));
    T2(i,1) = j;
    
end

% Rank and plot
%--------------------------------------------------------------------------
[s,j1] = sort(T1);
[s,j2] = sort(T2);

subplot(4,2,1), imagesc(1 - Y1(:,j1)')
title('reported cases: 1st peak','Fontsize',14), xlabel('days'), ylabel('local authority')
subplot(4,2,2), imagesc(1 - Y2(:,j2)')
title('reported cases: 2nd peak','Fontsize',14), xlabel('days')


% plot estimated reported cases
%--------------------------------------------------------------------------
n     = size(DC,1);
t1    = 1:200;
t2    = 201:n;
Y1    = zeros(numel(t1),numel(DCM));
Y2    = zeros(numel(t2),numel(DCM));
for i = 1:numel(DCM)
    
    % normalise to peak amplitude
    %----------------------------------------------------------------------
    Y       = DC;
    Y1(:,i) = Y(t1,i)/max(Y(t1,i));
    Y2(:,i) = Y(t2,i)/max(Y(t2,i));
    
    % amplitude and time of peaks
    %----------------------------------------------------------------------
    [m,j]   = max(Y(t1,i));
    T1(i,1) = j;
    [m,j]   = max(Y(t2,i));
    T2(i,1) = j;
    
end

% Rank and plot
%----------------------------------------------------------------------
[s,j1] = sort(T1);
[s,j2] = sort(T2);

subplot(4,2,3), imagesc(1 - Y1(:,j1)')
title('predicted reports: 1st peak','Fontsize',14), xlabel('days'), ylabel('local authority')
subplot(4,2,4), imagesc(1 - Y2(:,j2)')
title('predicted reports: 2nd peak','Fontsize',14), xlabel('days')


% plot estimated prevalence
%--------------------------------------------------------------------------
n     = size(DC,1);
t2    = 201:n;
Y2    = zeros(numel(t2),numel(DCM));
for i = 1:numel(DCM)
    
    % normalise to peak amplitude
    %----------------------------------------------------------------------
    Y       = DP;
    Y2(:,i) = Y(t2,i)/max(Y(t2,i));    
    [m,j]   = max(Y(t2,i));
    T2(i,1) = j;
    
end

% Rank and plot
%--------------------------------------------------------------------------
[s,j2] = sort(T2);
subplot(4,2,5), imagesc(1 - Y2(:,j2)')
title('predicted prevalence: 2nd peak','Fontsize',14), xlabel('days')

% plot estimated prevalence
%--------------------------------------------------------------------------
for i = 1:numel(DCM)
    
    % normalise to peak amplitude
    %----------------------------------------------------------------------
    Y       = DM;
    Y2(:,i) = Y(t2,i)/max(Y(t2,i));    
    [m,j]   = max(Y(t2,i));
    T2(i,1) = j;
end

% Rank and plot
%--------------------------------------------------------------------------
[s,j2] = sort(T2);
subplot(4,2,6), imagesc(1 - Y2(:,j2)')
title('no contact-rate fluctuations','Fontsize',14), xlabel('days')

% plot estimated prevalence
%--------------------------------------------------------------------------
for i = 1:numel(DCM)
    
    % normalise to peak amplitude
    %----------------------------------------------------------------------
    Y       = DV;
    Y2(:,i) = Y(t2,i)/max(Y(t2,i));    
    [m,j]   = max(Y(t2,i));
    T2(i,1) = j;
    
end

% Rank and plot
%--------------------------------------------------------------------------
[s,j2] = sort(T2);
subplot(4,2,7), imagesc(1 - Y2(:,j2)')
title('no transmission-risk fluctuations','Fontsize',14), xlabel('days')

% plot estimated prevalence
%--------------------------------------------------------------------------
n     = size(DC,1);
t2    = 201:n;
Y2    = zeros(numel(t2),numel(DCM));
for i = 1:numel(DCM)
    
    % normalise to peak amplitude
    %----------------------------------------------------------------------
    Y       = DB;
    Y2(:,i) = Y(t2,i)/max(Y(t2,i));    
    [m,j]   = max(Y(t2,i));
    T2(i,1) = j;
    
end

% Rank and plot
%--------------------------------------------------------------------------
[s,j2] = sort(T2);
subplot(4,2,8), imagesc(1 - Y2(:,j2)')
title('neither: second peak','Fontsize',14), xlabel('days')



%% regional analysis (two)
%==========================================================================

% get testing rates
%--------------------------------------------------------------------------
load DCM_UK
PCR = DCM.Y(:,4);
d0  = datenum(DCM.M.date,'dd-mm-yyyy');

% get regional data and parameter estimates
%--------------------------------------------------------------------------
load COVID_LA
PCR  = PCR/max(PCR);
dd   = dates(1) - d0;

% Rank regions by first and second peaks (before and after 200 days)
%==========================================================================
spm_figure('GetWin','Variation'), clf
%--------------------------------------------------------------------------
n     = numel(dates);
t1    = 1:200;
t2    = 201:n;
YC    = zeros(n,numel(DCM));
YD    = zeros(n,numel(DCM));
for i = 1:numel(DCM)
    
    % vectorised data and parameters
    %----------------------------------------------------------------------
    EP(:,i)  = spm_vec(DCM(i).Ep);
    
    % deal with missing data
    %----------------------------------------------------------------------
    S        = DCM(i).S;
    S(isnan(S)) = 0;
    j        = (1:size(S,1)) + n - size(S,1);
    YC(j,i)  = S(:,1);                     % new cases
    YD(j,i)  = S(:,2);                     % new deaths
    
end


% plot new cases, deaths and estimated prevalence
%--------------------------------------------------------------------------
tstr  = {'reported cases','reported deaths','estimated cases'};
xY    = {YC,YD,DP};
for k = 1:numel(xY)
    for i = 1:numel(DCM)
        
        % normalise to peak amplitude
        %------------------------------------------------------------------
        Y       = xY{k};
        Y1(:,i) = Y(t1,i)/max(Y(t1,i));
        Y2(:,i) = Y(t2,i)/max(Y(t2,i));
        
        % amplitude and time of peaks
        %------------------------------------------------------------------
        [m,j]   = max(Y(t1,i));
        M1(i,1) = m;
        T1(i,1) = j;
        
        [m,j]   = max(Y(t2,i));
        M2(i,1) = m;
        T2(i,1) = j;
        
    end

    % Rank and plot
    %----------------------------------------------------------------------
    [s,j1] = sort(T1);
    [s,j2] = sort(T2);
    
    subplot(4,2,(2*k - 1)), imagesc(1 - Y1(:,j1)')
    title(tstr{k},'Fontsize',14), xlabel('days'), ylabel('local authority')
    subplot(4,2,2*k),       imagesc(1 - Y2(:,j2)')
    title('second peak',   'Fontsize',14), xlabel('days')

end

% add testing
%--------------------------------------------------------------------------
subplot(4,2,1); hold on
plot((numel(DCM)  - PCR(t1)*300),'g')
subplot(4,2,2); hold on
plot((numel(DCM)  - PCR(t2)*300),'g')

% canonical correlation analysis
%==========================================================================
ip    = find(spm_vec(DCM(1).M.pC));
ip    = ip(ip <= 14);
param = fieldnames(DCM(1).M.pC);
param = param(ip);

% CCA
%-------------------------------------------------------------------------
Y     = [M1 T1 M2 T2];
c     = eye(numel(ip),numel(ip));
X     = EP(ip,:)';
X0    = ones(numel(DCM),1);
CVA   = spm_cva(Y,X,X0)

subplot(4,2,7), hold off
for i = 1:1
    plot(CVA.v(:,i),CVA.w(:,i),'.'), hold on
end
title('canonical variates','Fontsize',14), xlabel('parameters'), ylabel('timing')

subplot(4,2,8),   bar(CVA.W(:,1))
title('canonical vector','Fontsize',14),   xlabel('parameter'), ylabel('weight')
set(gca,'XTicklabel',param)

% correlations among amplitude and timing of first and second peaks
%--------------------------------------------------------------------------
[rho,pval] = corr([M1 T1 M2 T2])













