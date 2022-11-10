function [DCM] = DEM_COVID_DASH
% FORMAT [DCM] = DEM_COVID_DASH
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% The estimates of the reproduction rate and associated prevalence of
% infection in each region are based on a dynamic causal model of the
% coronavirus outbreak. This kind of modelling is distinct from
% conventional epidemiological modelling because the dynamic causal model
% incorporates things like self-isolation, social distancing, the
% probability of being tested and waiting for test results. This allows us
% to use regional reports of COVID-19 related deaths and new cases to model
% regional outbreaks.
% 
% In brief, the model assumes each region experiences an epidemic with
% similar characteristics but with different parameters, such as the number
% of contacts at home or work. And different mixtures of people who are
% more or less likely to catch (or transmit) the virus. These parameters
% are estimated from regional data and are then used to nowcast and
% forecast the evolution of the outbreak in terms of underlying or latent
% causes (such as the prevalence of infection). The effective population
% size is the number of people who are caught up in the outbreak, some of
% whom will be resistant to catching the virus. Of those that are not, some
% will become contagious and infect other people. From these estimates it
% is possible to evaluate the effective reproduction ratio at any point in
% time during the course of the outbreak, in addition to other quantitative
% estimates, such as the number of people currently infected or new cases
% of infection every day (that may or may not be identified).
% 
% The ensuing predictions complement equivalent estimates from
% epidemiological modelling based upon the history of outcomes observed so
% far. See https://www.mrc-bsu.cam.ac.uk/now-casting/ for a
% state-of-the-art transmission model. In principle, it is possible to
% compare the quality of dynamic causal and epidemiological models in terms
% of their model evidence or marginal likelihood. However, at the present
% time, it is difficult to estimate the evidence for epidemiological
% models; thereby precluding (Bayesian) model comparison.
% 
% Technical details about the dynamic causal model used here can be found
% at https://www.fil.ion.ucl.ac.uk/spm/covid-19/.
% 
% The (annotated) open source code creating these graphics is
% DEM_COVID_DASH.m

%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% get data
%==========================================================================
% N       = importdata('Region_UTLA_MetaData.csv');
% 
% regional populations (from Wikipedia)
%--------------------------------------------------------------------------
% South East 9,133,625 
% London  8,908,081 
% North West  5,479,093  
% East of England  6,201,214 
% West Midlands  5,900,757 
% South West  5,599,735 
% Yorkshire and the Humber  5,479,615  
% East Midlands  4,804,149   
% North East  2,657,909 
%--------------------------------------------------------------------------
Pop(1) = 6.201;          %     {'East Of England'         }
Pop(2) = 8.908;          %     {'London'                  }
Pop(3) = 5.900 + 4.804;  %     {'Midlands'                }.
Pop(4) = 7.292 + 5.479;  %     {'North East And Yorkshire'}
Pop(5) = 5.479;          %     {'North West'              }
Pop(6) = 9.133;          %     {'South East'              }
Pop(7) = 5.599;          %     {'South West'              }

% retrieve recent data from https://coronavirus.data.gov.uk
%--------------------------------------------------------------------------
url  = 'https://coronavirus.data.gov.uk/downloads/csv/';
websave('coronavirus-cases_latest.csv',[url,'coronavirus-cases_latest.csv']);

% retrieve recent data from https://www.england.nhs.uk
%--------------------------------------------------------------------------
URL   = 'https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/';
URL   = [URL datestr(date,'yyyy/mm') '/'];
for i = 0:4
    try
        dstr = datestr(datenum(date) - i,'dd-mmmm-yyyy');
        if strcmp(dstr(1),'0'),dstr = dstr(2:end); end
        url  = [URL 'COVID-19-total-announced-deaths-' dstr '.xlsx'];
        fprintf('Trying %s\n',url);
        websave('COVID-19-total-announced-deaths.xlsx',url);
        fprintf('Using %s\n',url);
        break
    end
end

% load data
%--------------------------------------------------------------------------
C    = importdata('coronavirus-cases_latest.csv');
D    = importdata('COVID-19-total-announced-deaths.xlsx');
L    = importdata('Local_Authority_District_to_Region_(April_2019)_Lookup_in_England.xlsx');

DN   = D.textdata.Tab1DeathsByRegion(15,4:end - 4);
DR   = D.textdata.Tab1DeathsByRegion(18:end,1);

DY   = D.data.Tab1DeathsByRegion(3:end,2:end - 4);
DY   = D.data.Tab2Deaths0x2DNoPostTest(3:end,2:end - 4) + DY;

% find death by date
%--------------------------------------------------------------------------
for i = 1:numel(DN)
    try
        dn(i)  = datenum(DN(i),'dd-mmm-yy');
    catch
        break;
    end
end
i       = find(dn);
dn      = dn(i);
DY      = DY(:,i);

% match death rates with reporting new cases
%--------------------------------------------------------------------------
regions = unique(L.textdata(2:end,5));

% get cumulative cases for each region
%--------------------------------------------------------------------------
for r = 1:numel(regions)
    i     = find(ismember(L.textdata(:,5),regions{r}));
    i     = find(ismember(C.textdata(:,2),L.textdata(i,2)));
    
    dat   = unique(C.textdata(i,4));
    cases = C.data(i - 1,3);
    for d = 1:numel(dat)
       j        = find(ismember(C.textdata(i,4),dat(d)));
       cy{r}(d) = sum(cases(j));
       cn{r}(d) = datenum(dat(d),'yyyy-mm-dd');
    end
end

% and associate with death rates
%--------------------------------------------------------------------------
for i = 1:numel(dn)
    for r = 1:numel(regions)
        d       = abs(cn{r} - dn(i));
        [d,j]   = min(d);
        CY(r,i) = cy{r}(j);
    end
end

% combine regional data into NHS regions
%--------------------------------------------------------------------------
nr    = numel(DR);
for r = 1:nr
    if strcmp(DR{r},'Midlands')
        j = [1,8];
    elseif strcmp(DR{r},'North East And Yorkshire')
        j = [4,9];
    elseif strcmp(DR{r},'East Of England')
        j = 2;
    else
        j = find(ismember(regions,DR(r)));
    end
    CR(r,:) = sum(CY(j,:),1);
end

% prepend 16 days to timeseries
%--------------------------------------------------------------------------
DY   = [zeros(nr,16), DY];
CR   = [zeros(nr,16), CR];
dn   = [(dn(1) - flip(1:16)) dn];


% fit each regional dataset
%==========================================================================
for r = 1:numel(DR)
    
    % Get data for this region
    %======================================================================
    s        = 16;
    Y        = [spm_hist_smooth(DY(r,:),s), ...
                spm_hist_smooth(gradient(CR(r,:)),s)];
    
    % get (Gaussian) priors over model parameters
    %----------------------------------------------------------------------
    [pE,pC]  = spm_SARS_priors;
    
    % priors for this analysis
    %----------------------------------------------------------------------
    pE.N   = log(Pop(r));                 % population of region (M)
    pC.N   = 0;
    pE.n   = 4;                           % initial number of cases (n)

    
    % variational Laplace (estimating log evidence (F) and posteriors)
    %======================================================================
    
    % complete model specification
    %----------------------------------------------------------------------
    M.date = datestr(dn(1),'dd-mm-yyyy'); % date of first time point
    M.G    = @spm_SARS_gen;               % generative function
    M.FS   = @(Y)sqrt(Y);                 % feature selection  (link function)
    M.pE   = pE;                          % prior expectations (parameters)
    M.pC   = pC;                          % prior covariances  (parameters)
    M.hE   = 2;                           % prior expectation  (log-precision)
    M.hC   = 1/512;                       % prior covariances  (log-precision)
    M.T    = size(Y,1);                   % number of samples
    U      = [1 2];                       % outputs to model
    
    % model inversion with Variational Laplace (Gauss Newton)
    %----------------------------------------------------------------------
    [Ep,Cp]   = spm_nlsi_GN(M,U,Y);
    
    % save prior and posterior estimates (and log evidence)
    %----------------------------------------------------------------------
    DCM(r).M  = M;
    DCM(r).Ep = Ep;
    DCM(r).Cp = Cp;
    DCM(r).Y  = Y;
    
    % now-casting for this region and date
    %======================================================================
    STR = [DR{r} '   ' datestr(dn(end))];
    
    % show predictions with confidence intervals
    %----------------------------------------------------------------------
    spm_figure('GetWin',STR); clf;
    %----------------------------------------------------------------------
    M.T = 385;
    t   = datenum(date);
    spm_SARS_ci(DCM(r).Ep,DCM(r).Cp,DCM(r).Y,[1 2 4],M);
    
    subplot(3,2,5), XLim = get(gca,'XLim'); YLim = get(gca,'YLim');
    hold on, plot(XLim,[1,1],'r-.'), plot([t,t],YLim,'b:')
    title('Reproduction ratio (R)','Fontsize',16)
    ylabel('reproduction ratio')
    
    % supplement with table of posterior expectations
    %----------------------------------------------------------------------
    subplot(2,2,2), title('Cumulative deaths','Fontsize',16), ylabel('deaths')
    subplot(3,2,1), ylabel('deaths per day')
    %----------------------------------------------------------------------
    % Y(:,1) - number of new deaths
    % Y(:,2) - number of new cases
    % Y(:,3) - CCU bed occupancy
    % Y(:,4) - effective reproduction rate (R)
    % Y(:,5) - herd immunity
    % Y(:,6) - total number of tests
    % Y(:,7) - contagion risk (%)
    % Y(:,8) - prevalence of infection (%)
    % Y(:,9) - number of infected at home, untested and asymptomatic
    % and plot latent or hidden states
    %----------------------------------------------------------------------
    [R,X] = spm_SARS_gen(DCM(r).Ep,DCM(r).M,[4 5 8 9]);
    
    subplot(2,2,4), cla reset, axis([0 1 0 1])
    title(STR,'FontSize',16)
    
    Tab(r,1) = Pop(r);
    str      = sprintf('Census population %.2f million',Tab(r,1));
    text(0,1.0,str,'FontSize',14,'Color','k')
    
    Tab(r,2) = exp(Ep.N)*exp(Ep.o);
    str      = sprintf('Initially exposed %.2f million',Tab(r,2));
    text(0,0.9,str,'FontSize',14,'Color','k')
    
    Tab(r,3) = R(end,1);
    str      = sprintf('Reproduction ratio %.2f',Tab(r,3));
    text(0,0.8,str,'FontSize',14,'Color','k')
    
    Tab(r,4) = R(end,4);
    str      = sprintf('Infected, asymptomatic people %.0f',Tab(r,4));
    text(0,0.7,str,'FontSize',14,'Color','k')
    
    Tab(r,5) = 1e4*exp(Ep.N)*R(end,3)/exp(Ep.Tin);
    str      = sprintf('New infections today %.0f',Tab(r,5));
    text(0,0.6,str,'FontSize',14,'Color','k')
    
    Tab(r,6) = R(end,3)*exp(Ep.N)/Pop(r);
    str      = sprintf('Prevalence of infection %.2f%s',Tab(r,6),'%');
    text(0,0.5,str,'FontSize',14,'Color','k')
    
    Tab(r,7) = R(end,2)*exp(Ep.N)/Pop(r);
    str = sprintf('Prevalence of immunity %.1f%s',Tab(r,7),'%');
    text(0,0.4,str,'FontSize',14,'Color','k')
    
    str = {'The prevalences refer to the total population' ...
           'as estimated from the new cases and deaths (shown as' ...
           'black dots on the upper left panels)'};
    text(0,0.2,str,'FontSize',10,'Color','k')
    
    str = {'The new infections today can be regarded as the minimum' ...
           'number of people who need to be found on a daily basis' ...
           'by an FTTIS program'};
    text(0,0.0,str,'FontSize',10,'Color','m')
    
    spm_axis tight, axis off
    
end

% summary page with reproduction rate and prevalence of infection
%--------------------------------------------------------------------------
spm_figure('GetWin','Summary'); clf;

subplot(2,2,3), bar(Tab(:,3))
set(gca,'XTickLabel',DR,'XTickLabelRotation',90)
ylabel('Reproduction rate')
title('Reproduction rate (R)','Fontsize',16)
set(gca,'YLim',[.5 1]); axis square, box off

subplot(2,2,1), bar(Tab(:,6))
ylabel('Percent of population')
set(gca,'XTickLabel',DR,'XTickLabelRotation',90)
title('Prevalence of infection','Fontsize',16)
axis square, box off

% notes
%--------------------------------------------------------------------------
subplot(2,2,2),cla, axis([0 1 0 1])
h  = help('DEM_COVID_DASH');
text(-.2,1.2,h,'FontSize',9,'Color','k','VerticalAlignment','cap')
spm_axis tight, axis off





