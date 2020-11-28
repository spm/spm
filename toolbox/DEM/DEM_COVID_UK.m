function DCM = DEM_COVID_UK
% FORMAT DCM = DEM_COVID_UK
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% This routine illustrates Bayesian model comparison using a line search
% over periods of imunity and pooling over countries. In brief,32 countries
% are inverted and 16 with the most informative posterior over the period
% of immunity are retained for Bayesian parameter averaging. The Christian
% predictive densities are then provided in various formats for the average
% country and (16) individual countries.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DEM_COVID_UK.m 8024 2020-11-28 12:09:53Z karl $

% set up and preliminaries
%==========================================================================
% https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
% https://www.ndm.ox.ac.uk/covid-19/covid-19-infection-survey/results
% https://coronavirus.data.gov.uk/
% https://covid.joinzoe.com/data#levels-over-time
% https://www.gov.uk/guidance/the-r-number-in-the-uk#history
% https://www.gov.uk/government/statistics/transport-use-during-the-coronavirus-covid-19-pandemic
% https://www.google.com/covid19/mobility/

% get figure and data
%--------------------------------------------------------------------------
Fsi       = spm_figure('GetWin','SI'); clf;

cd('C:\Users\karl\Dropbox\Coronavirus\Dashboard')

url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newCasesByPublishDate&format=csv';
writetable(webread(url),'cases.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newDeaths28DaysByDeathDate&format=csv';
writetable(webread(url),'deaths.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=covidOccupiedMVBeds&format=csv';
writetable(webread(url),'critical.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newPillarOneTwoTestsByPublishDate&format=csv';
writetable(webread(url),'tests.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newAdmissions&format=csv';
writetable(webread(url),'admissions.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newOnsDeathsByRegistrationDate&format=csv';
writetable(webread(url),'certified.csv');

cases      = importdata('cases.csv');
deaths     = importdata('deaths.csv');
ccu        = importdata('critical.csv');
tests      = importdata('tests.csv');
certified  = importdata('certified.csv');
admissions = importdata('admissions.csv');

serology   = importdata('seropositive.csv');
survey     = importdata('survey.csv');
symptoms   = importdata('symptoms.csv');
ratio      = importdata('rate.csv');
mobility   = importdata('mobility.csv');
transport  = importdata('transport.csv');


dstr      = 'dd/mm/yyyy';
dONS      = 'yyyy-mm-dd';

% created data structure
%--------------------------------------------------------------------------
Y(1).type = 'PCR cases (ONS)'; % daily PCR positive cases (by specimen)
Y(1).unit = 'number/day';
Y(1).U    = 2;
Y(1).date = datenum(cases.textdata(2:end,1),dONS);
Y(1).Y    = cases.data(:,1);
Y(1).h    = 0;

Y(2).type = 'Daily deaths (ONS: 28-days)'; % daily covid-related deaths (28 days)
Y(2).unit = 'number/day';
Y(2).U    = 1;
Y(2).date = datenum(deaths.textdata(2:end,1),dONS);
Y(2).Y    = deaths.data(:,1);
Y(2).h    = 2;

Y(3).type = 'Ventilated patients (ONS)'; % CCU occupancy (mechanical)
Y(3).unit = 'number';
Y(3).U    = 3;
Y(3).date = datenum(ccu.textdata(2:end,1),dONS);
Y(3).Y    = ccu.data(:,1);
Y(3).h    = 0;

Y(4).type = 'PCR tests (ONS)'; % daily PCR tests performed
Y(4).unit = 'number/day';
Y(4).U    = 6;
Y(4).date = datenum(tests.textdata(2:end,1),dONS);
Y(4).Y    = tests.data(:,1);
Y(4).h    = 0;

Y(5).type = 'Prevalence (ONS)'; % number of people infected (England)
Y(5).unit = 'number';
Y(5).U    = 11;
Y(5).date = datenum(survey.textdata(2:end,1),dstr) - 7;
Y(5).Y    = survey.data(:,1)*66.79/56.28;
Y(5).h    = 0;

Y(6).type = 'Seropositive (GOV)'; % percentage seropositive
Y(6).unit = 'percent';
Y(6).U    = 5;
Y(6).date = datenum(serology.textdata(2:end,1),dstr);
Y(6).Y    = serology.data(:,1);
Y(6).h    = 0;

Y(7).type = 'Symptoms (KCL)'; % number of people reporting symptoms (UK)
Y(7).unit = 'number';
Y(7).U    = 12;
Y(7).date = datenum(symptoms.textdata(2:end,1),dstr);
Y(7).Y    = symptoms.data(:,1);
Y(7).h    = 0;

Y(8).type = 'R-ratio (GOV)'; % the production ratio
Y(8).unit = 'ratio';
Y(8).U    = 4;
Y(8).date = [datenum(ratio.textdata(2:end,1),dstr) - 13; ...
             datenum(ratio.textdata(2:end,1),dstr) - 14];
Y(8).Y    = [ratio.data(:,1); ratio.data(:,2)];
Y(8).h    = 0;

Y(9).type = 'Transport (GOV)'; % cars (percent)
Y(9).unit = 'percent';
Y(9).U    = 13;
Y(9).date = datenum(transport.textdata(2:end,1),dstr);
Y(9).Y    = transport.data(:,1)*100;
Y(9).h    = 0;

Y(10).type = 'Work (Google)'; % work (percent)
Y(10).unit = 'percent';
Y(10).U    = 14;
Y(10).date = datenum(mobility.textdata(2:end,1),dstr);
Y(10).Y    = mobility.data(:,5) + 100;
Y(10).h    = 0;

Y(11).type = 'Certified deaths (ONS)'; % weekly covid related deaths
Y(11).unit = 'number';
Y(11).U    = 15;
Y(11).date = datenum(certified.textdata(2:end,1),dONS) - 7;
Y(11).Y    = certified.data(:,1)/7;
Y(11).h    = 0;

Y(12).type = 'Admissions (ONS)'; % admissions to hospital
Y(12).unit = 'number';
Y(12).U    = 16;
Y(12).date = datenum(admissions.textdata(2:end,1),dONS);
Y(12).Y    = admissions.data(:,1);
Y(12).h    = 0;

% remove NANs, smooth and sort by date
%==========================================================================
M.date  = datestr(min(spm_vec(Y.date)),'dd-mm-yyyy');
[Y,YS]  = spm_COVID_Y(Y,M.date);

% data structure with vectorised data and covariance components
%--------------------------------------------------------------------------
xY.y    = spm_vec(Y.Y);
xY.Q    = spm_Ce([Y.n]);
hE      = spm_vec(Y.h);

% get and set priors
%--------------------------------------------------------------------------
[pE,pC] = spm_SARS_priors;
pE.N    = log(66.65);
pC.N    = 0;
pE.n    = -8;

% coefficients for likelihood model
%--------------------------------------------------------------------------
pE.mo   = log([2,1/2]);         % coefficients for mobility
pC.mo   = exp([-2,-2]);         % prior variance

pE.wo   = log([2,1/2]);         % coefficients for workplace
pC.wo   = exp([-2,-2]);         % prior variance

pE.ho   = log([1,  1]);         % coefficients for admissions
pC.ho   = exp([-2,-2]);         % prior variance

pE.dc   = log(1);               % death (28 days) coefficients
pC.dc   = exp(-2);              % prior variance

% model specification
%==========================================================================
M.Nmax  = 64;                   % maximum number of iterations
M.G     = @spm_SARS_gen;        % generative function
M.FS    = @(Y)real(sqrt(Y));    % feature selection  (link function)
M.pE    = pE;                   % prior expectations (parameters)
M.pC    = pC;                   % prior covariances  (parameters)
M.hE    = hE;                   % prior expectation  (log-precision)
M.hC    = 1/512;                % prior covariances  (log-precision)
M.T     = Y;                    % data structure
U       = spm_vec(Y.U);         % outputs to model

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
[Ep,Cp,Eh] = spm_nlsi_GN(M,U,xY);
DCM.M   = M;
DCM.Ep  = Ep;
DCM.Eh  = Eh;
DCM.Cp  = Cp;
DCM.Y   = YS;
DCM.U   = U;


% posterior predictions
%==========================================================================
spm_figure('GetWin','United Kingdom'); clf;
%--------------------------------------------------------------------------
M.T     = datenum('01-04-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
U       = U(1:min(numel(U),3));
[Z,X]   = spm_SARS_gen(DCM.Ep,M,U);
spm_SARS_plot(Z,X,YS,[],U)
U       = DCM.U;

spm_figure('GetWin','outcomes (1)'); clf;
%--------------------------------------------------------------------------
j     =  0;
for i = 1:numel(Y)
    j = j + 1;
    subplot(4,2,j)
    spm_SARS_ci(Ep,Cp,YS(:,i),U(i),M);
    title(Y(i).type,'FontSize',14), ylabel(Y(i).unit)
    
    % add R = 1 and current date
    %----------------------------------------------------------------------
    if i == 8
        plot(get(gca,'XLim'),[1,1],'-.r')
        plot(datenum(date)*[1,1],get(gca,'YLim'),'-.b')
        g  = get(gca,'YLim');
        plot(datenum(date)*[1,1],g,'-.b')
        set(gca,'YLim',[0 min(4,g(2))]), ylabel('ratio')
    end
    
    % new figure
    %----------------------------------------------------------------------
    if j > 7
        spm_figure('GetWin','outcomes (2)'); clf;
        j = 0;
    end
end

% reproduction ratio
%--------------------------------------------------------------------------
if numel(Y) < 8
    j = j + 1;
    subplot(4,2,j)
    spm_SARS_ci(Ep,Cp,[],4,M);
    plot(get(gca,'XLim'),[1,1],'-.r')
    g  = get(gca,'YLim');
    plot(datenum(date)*[1,1],g,'-.b')
    set(gca,'YLim',[0 min(4,g(2))]), ylabel('ratio')
end

% save figures
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes (2)');
savefig(gcf,'Fig2')

spm_figure('GetWin','outcomes (1)');
savefig(gcf,'Fig1')

spm_figure('GetWin','United Kingdom');
savefig(gcf,'Fig3')


% Table
%--------------------------------------------------------------------------
Tab = spm_COVID_table(Ep,Cp,M)

cd('C:\Users\karl\Dropbox\Coronavirus')

return

% repeat with rapid loss of immunity
%==========================================================================

% spm_figure('GetWin',[country ': 6 months']); clf;
% %------------------------------------------------------------------------
% Ep.Tim   = log(6);
% [Z,X]    = spm_SARS_gen(Ep,M,U(1:3));
% spm_SARS_plot(Z,X,YS)
% 
% spm_figure('GetWin','death rates'); hold on
% %------------------------------------------------------------------------
% spm_SARS_ci(Ep,Cp,[],1,M);
% Ep.Tim = DCM.Ep.Tim;

% repeat with efficient FTTIS
%==========================================================================
spm_figure('GetWin','death rates'); hold on
M.FTT = 1/4;
for i = [0 4 8]*7
    M.TTT = datenum(date) - datenum(M.date,'dd-mm-yyyy') + i;
    spm_SARS_ci(Ep,Cp,[],1,M);
end

% repeat for hospital admissions
%==========================================================================
spm_figure('GetWin','CCU admissions'); hold on
%--------------------------------------------------------------------------
for i = [0 4 8]*7
    M.TTT = datenum(date) - datenum(M.date,'dd-mm-yyyy') + i;
    spm_SARS_ci(Ep,Cp,[],3,M);
end

% repeat for Circuit break
%==========================================================================
spm_figure('GetWin','Circuit break'); clf
%--------------------------------------------------------------------------
CBT   = datenum(date) - datenum(DCM.M.date,'dd-mm-yyyy') + 4;
CBD   = 14;                         % duration of circuit breaker

M     = DCM.M;
M.T   = datenum('01-1-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
[Z,X] = spm_SARS_gen(DCM.Ep,M,[2 1 3]);
spm_SARS_plot(Z,X,DCM.YS,[],[2 1 3])
for j = 1:6
    subplot(3,2,j), hold on
    set(gca,'ColorOrderIndex',1);
end

M.CBT = CBT;
M.CBD = CBD;
[Z,X] = spm_SARS_gen(DCM.Ep,M,[2 1 3]);
spm_SARS_plot(Z,X,DCM.YS)

% fatalities in confidence intervals
%--------------------------------------------------------------------------
spm_figure('GetWin','Circuit break - deaths'); clf
M     = DCM.M;
M.T   = datenum('01-01-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
spm_SARS_ci(DCM.Ep,DCM.Cp,DCM.YS,1,M); hold on

M.CBT = CBT;
M.CBD = CBD;
spm_SARS_ci(DCM.Ep,DCM.Cp,[],1,M);


M.CBT = CBT;
M.CBD = CBD;
M.FTT = 1/4;
M.TTT = CBT + 28;
spm_SARS_ci(DCM.Ep,DCM.Cp,[],1,M);

M.CBT = CBT;
M.CBD = CBD;
M.FTT = .8;
M.TTT = CBT + 28;
spm_SARS_ci(DCM.Ep,DCM.Cp,[],1,M);

t   = datenum('01-1-2021','dd-mm-yyyy');
t0  = t - 7*24;
set(gca,'XLim',[t0,t])
set(gca,'YLim',[0,200])
set(gca,'XTick',[t0:14:t])
datetick('x','mmm-dd','keeplimits','keepticks')




