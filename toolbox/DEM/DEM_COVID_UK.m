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
% $Id: DEM_COVID_UK.m 8036 2020-12-20 19:19:56Z karl $

% set up and preliminaries
%==========================================================================
% https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
% https://www.ndm.ox.ac.uk/covid-19/covid-19-infection-survey/results
% https://coronavirus.data.gov.uk/
% https://covid.joinzoe.com/data#levels-over-time
% https://www.gov.uk/guidance/the-r-number-in-the-uk#history
% https://www.gov.uk/government/statistics/transport-use-during-the-coronavirus-covid-19-pandemic
% https://www.google.com/covid19/mobility/

% Age group		Yes	No	Unkown presence of pre-existing condition 	Total
% ________________________________________
% Total		    37,847	1,723	0	39,570	
% 0 - 19 		19	4	0	23	
% 20 - 39		228	42	0	270	
% 40 - 59		2,527	303	0	2,830	
% 60 - 79		14,417	713	0	15,130	
% 80+		    20,656	661	0	21,317	



% web options
%--------------------------------------------------------------------------
options = weboptions('ContentType','table'); 
options.Timeout = 20;

% England to UK operation conversion
%--------------------------------------------------------------------------
EnglandUK     = 66.79/56.28;
EngandWalesUK = 66.79/(56.28 + 3.15);

% get figure and data
%--------------------------------------------------------------------------
spm_figure('GetWin','SI'); clf;

cd('C:\Users\karl\Dropbox\Coronavirus\Dashboard')

url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newCasesBySpecimenDate&format=csv';
writetable(webread(url,options),'cases.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newDeaths28DaysByDeathDate&format=csv';
writetable(webread(url,options),'deaths.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=covidOccupiedMVBeds&format=csv';
writetable(webread(url,options),'critical.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newPillarOneTwoTestsByPublishDate&format=csv';
writetable(webread(url,options),'tests.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newAdmissions&format=csv';
writetable(webread(url,options),'admissions.csv');
url        = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newOnsDeathsByRegistrationDate&format=csv';
writetable(webread(url,options),'certified.csv');
url        = 'https://www.ons.gov.uk/generator?uri=/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/weekending13november2020/00ba3836&format=csv';
url        = 'https://www.ons.gov.uk/generator?uri=/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/weekending27november2020/96f0e889&format=csv';
writetable(webread(url,options),'place.csv');
url        = 'https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2020/11/COVID-19-total-announced-deaths-27-November-2020.xlsx';
url        = 'https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2020/12/COVID-19-total-announced-deaths-18-December-2020.xlsx';
[num,txt]  = xlsread(websave('ages.xlsx',url),5,'E16:KC23');


cases      = importdata('cases.csv');
deaths     = importdata('deaths.csv');
ccu        = importdata('critical.csv');
tests      = importdata('tests.csv');
certified  = importdata('certified.csv');
admissions = importdata('admissions.csv');

serology   = importdata('seropositive.csv');
survey     = importdata('survey.csv');
symptoms   = importdata('symptoms.csv');
ratio      = importdata('ratio.csv');
mobility   = importdata('mobility.csv');
transport  = importdata('transport.csv');
place      = importdata('place.csv');

% created data structure
%--------------------------------------------------------------------------
Y(1).type = 'PCR cases (ONS)'; % daily PCR positive cases (by specimen)
Y(1).unit = 'number/day';
Y(1).U    = 2;
Y(1).date = datenum(cases.textdata(2:end,1),'yyyy-mm-dd');
Y(1).Y    = cases.data(:,1);
Y(1).h    = 0;

Y(2).type = 'Daily deaths (ONS: 28-days)'; % daily covid-related deaths (28 days)
Y(2).unit = 'number/day';
Y(2).U    = 1;
Y(2).date = datenum(deaths.textdata(2:end,1),'yyyy-mm-dd');
Y(2).Y    = deaths.data(:,1);
Y(2).h    = 0;

Y(3).type = 'Ventilated patients (ONS)'; % CCU occupancy (mechanical)
Y(3).unit = 'number';
Y(3).U    = 3;
Y(3).date = datenum(ccu.textdata(2:end,1),'yyyy-mm-dd');
Y(3).Y    = ccu.data(:,1);
Y(3).h    = 0;

Y(4).type = 'PCR tests (ONS)'; % daily PCR tests performed
Y(4).unit = 'number/day';
Y(4).U    = 6;
Y(4).date = datenum(tests.textdata(2:end,1),'yyyy-mm-dd');
Y(4).Y    = tests.data(:,1);
Y(4).h    = 0;

Y(5).type = 'Prevalence (ONS)'; % number of people infected (England)
Y(5).unit = 'number';
Y(5).U    = 11;
Y(5).date = datenum(survey.textdata(2:end,1),'dd/mm/yyyy') - 7;
Y(5).Y    = survey.data(:,1)*EnglandUK;
Y(5).h    = 0;

Y(6).type = 'Seropositive (GOV)'; % percentage seropositive
Y(6).unit = 'percent';
Y(6).U    = 5;
Y(6).date = datenum(serology.textdata(2:end,1),'dd/mm/yyyy');
Y(6).Y    = serology.data(:,1);
Y(6).h    = 0;

Y(7).type = 'Symptoms (KCL)'; % number of people reporting symptoms (UK)
Y(7).unit = 'number';
Y(7).U    = 12;
Y(7).date = datenum(symptoms.textdata(2:end,1),'dd/mm/yyyy');
Y(7).Y    = symptoms.data(:,1);
Y(7).h    = 0;

Y(8).type = 'R-ratio (MRC/GOV)'; % the production ratio
Y(8).unit = 'ratio';
Y(8).U    = 4;
Y(8).date = [datenum(ratio.textdata(2:end,1),'dd/mm/yyyy') - 13; ...
             datenum(ratio.textdata(2:end,1),'dd/mm/yyyy') - 14];
Y(8).Y    = [ratio.data(:,1); ratio.data(:,2)];
Y(8).h    = 0;

Y(9).type = 'Transport (GOV)'; % cars (percent)
Y(9).unit = 'percent';
Y(9).U    = 13;
Y(9).date = datenum(transport.textdata(2:end,1),'dd/mm/yyyy');
Y(9).Y    = transport.data(:,1)*100;
Y(9).h    = 0;

Y(10).type = 'Retail (Google)'; % retail and recreation (percent)
Y(10).unit = 'percent';
Y(10).U    = 14;
Y(10).date = datenum(mobility.textdata(2:end,1),'dd/mm/yyyy');
Y(10).Y    = mobility.data(:,1) + 100;
Y(10).h    = 0;

Y(11).type = 'Certified deaths (ONS)'; % weekly covid related deaths
Y(11).unit = 'number';
Y(11).U    = 15;
Y(11).date = datenum(certified.textdata(2:end,1),'yyyy-mm-dd') - 7;
Y(11).Y    = certified.data(:,1)/7;
Y(11).h    = 0;

Y(12).type = 'Admissions (ONS)'; % admissions to hospital
Y(12).unit = 'number';
Y(12).U    = 16;
Y(12).date = datenum(admissions.textdata(2:end,1),'yyyy-mm-dd');
Y(12).Y    = admissions.data(:,1);
Y(12).h    = 0;

Y(13).type = 'Hospital deaths (PHE)'; % hospital deaths
Y(13).unit = 'number';
Y(13).U    = 17;
Y(13).date = datenum(place.textdata(2:end,1),'dd mmm');
Y(13).Y    = place.data(:,4)*EngandWalesUK;
Y(13).h    = 0;

Y(14).type = 'Hospital/Other deaths (PHE)'; % nonhospital deaths
Y(14).unit = 'number';
Y(14).U    = 18;
Y(14).date = datenum(place.textdata(2:end,1),'dd mmm') - 7;
Y(14).Y    = sum(place.data(:,1:3),2)*EngandWalesUK;
Y(14).h    = 0;

Y(15).type = 'deaths > 60 (PHE)'; % deaths (English hospitals)
Y(15).unit = 'number';
Y(15).U    = 19;
Y(15).date = datenum(txt,'dd/mm/yyyy');
Y(15).Y    = sum(num(6:7,:))'*EnglandUK;
Y(15).h    = 0;

Y(16).type = 'deaths </> 60 (PHE)'; % deaths (English hospitals)
Y(16).unit = 'number';
Y(16).U    = 20;
Y(16).date = datenum(txt,'dd/mm/yyyy');
Y(16).Y    = sum(num(3:5,:))'*EnglandUK;
Y(16).h    = 0;

% normalise total deaths to ONS certified deaths
%--------------------------------------------------------------------------
d       = certified.data(:,1);
N       = sum(d(isfinite(d)));

d       = Y(13).Y;
n       = sum(d(isfinite(d)));
d       = Y(14).Y;
n       = sum(d(isfinite(d))) + n;
Y(13).Y = Y(13).Y*N/n;
Y(14).Y = Y(14).Y*N/n;

d       = Y(15).Y;
n       = sum(d(isfinite(d)));
d       = Y(16).Y;
n       = sum(d(isfinite(d))) + n;
Y(15).Y = Y(15).Y*N/n;
Y(16).Y = Y(16).Y*N/n;


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
pE.N    = log(67.886);
pC.N    = 0;
pE.n    = -5;

% coefficients for likelihood model
%--------------------------------------------------------------------------
pE.dc   = log([1 1]);           % coefficients for death (28 days)
pC.dc   = [1 1]/8;              % prior variance
pE.mv   = log([1 1]);           % coefficients for ventilation
pC.mv   = [1 1]/8;              % prior variance
pE.ho   = log([1 1]);           % coefficients for admissions
pC.ho   = [1 1]/8;              % prior variance

pE.mo   = log([1.5,0.3]);       % coefficients for mobility
pC.mo   = [1 1];                % prior variance
pE.wo   = log([1.5,0.3]);       % coefficients for workplace
pC.wo   = [1 1];                % prior variance

pE.ag   = zeros(2,3);           % coefficients for age-related deaths
pC.ag   = ones(size(pE.ag));    % prior variance

% coefficients for mixture model
%--------------------------------------------------------------------------
% pE.D(1).o   = pE.o   - 1/4;
% pE.D(1).Nin = pE.Nin + 1/4;
% pE.D(1).Nou = pE.Nou + 1/4;
% pE.D(1).trn = pE.trn + 1/4;
% pE.D(1).trm = pE.trm + 1/4;
% 
% for i = 1:numel(pE.D)
%     param = fieldnames(pE.D(i));
%     for j = 1:numel(param)
%         pC.D(i).(param{j}) = pC.(param{j});
%     end
% end


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
%==========================================================================
[Ep,Cp,Eh] = spm_nlsi_GN(M,U,xY);


% fluctuations (adiabatic mean field approximation)
%--------------------------------------------------------------------------
fluct = {'vir','pcr'};
for f = 1:numel(fluct)
    
    i               = 1:size(Cp,1);          % number of parameters
    C               = Cp;                    % empirical prior covariance
    M.pE            = Ep;                    % empirical prior expectation
    M.pC            = spm_zeros(M.pC);       % fix expectations
    M.pE.(fluct{f}) = zeros(1,16);           % add new prior expectation
    M.pC.(fluct{f}) =  ones(1,16);           % add new prior covariance
    
    [Ep,Cp]  = spm_nlsi_GN(M,U,xY);          % new posterior expectation
    Cp(i,i)  = C;                            % new posterior covariance
    
end

% save in DCM structure
%--------------------------------------------------------------------------
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
M.T       = datenum('01-04-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
U         = U(1:min(numel(U),3));
[H,X,Z,R] = spm_SARS_gen(Ep,M,U);
spm_SARS_plot(H,X,YS,[],U)
U         = DCM.U;

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
        set(gca,'YLim',[0 min(5,g(2))]), ylabel('ratio')
    end
       
    % hold plot
    %----------------------------------------------------------------------
    if i == 13
        j = j - 1; hold on
    end
    
    % hold plot
    %----------------------------------------------------------------------
    if i == 15 
        j = j - 1; hold on
    end

    % new figure
    %----------------------------------------------------------------------
    if i == 8
        spm_figure('GetWin','outcomes (2)');clf
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
    set(gca,'YLim',[0 min(5,g(2))]), ylabel('ratio')
end

% infection fatality ratios (%)
%--------------------------------------------------------------------------
j = j + 1;
subplot(4,2,j)
spm_SARS_ci(Ep,Cp,[],21,M)
ylabel('percent'),  title('Infection fatality ratio','FontSize',14)

% transmission strength
%--------------------------------------------------------------------------
j    = j + 1;
subplot(4,2,j)
plot([R.Ptrn]), spm_axis tight
title('Transmission strength','FontSize',14)
set(gca,'YLim',[0 1]), ylabel('probability')
hold on, plot([1,1]*size(DCM.Y,1),[0,1],':'), hold off

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

save('DCM_UK.mat','DCM')
cd('C:\Users\karl\Dropbox\Coronavirus')
save('DCM_UK.mat','DCM')

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




