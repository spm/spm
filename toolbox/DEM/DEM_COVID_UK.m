function DCM = DEM_COVID_UK(fluct)
% FORMAT DCM = DEM_COVID_UK(fluct)
% fluct - fluctuations; e.g., fluct = {'mob','pcr'}
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
% $Id: DEM_COVID_UK.m 8039 2021-01-03 09:46:59Z karl $

% set up and preliminaries
%==========================================================================
% https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
% https://www.ndm.ox.ac.uk/covid-19/covid-19-infection-survey/results
% https://coronavirus.data.gov.uk/
% https://covid.joinzoe.com/data#levels-over-time
% https://www.gov.uk/guidance/the-r-number-in-the-uk#history
% https://www.gov.uk/government/statistics/transport-use-during-the-coronavirus-covid-19-pandemic
% https://www.google.com/covid19/mobility/

% web options
%--------------------------------------------------------------------------
if nargin < 1, fluct = {'vir','mob','pcr'}; end

% web options
%--------------------------------------------------------------------------
options = weboptions('ContentType','table');
options.Timeout = 20;

% England to UK operation conversion
%--------------------------------------------------------------------------
EnglandUK     = 66.79/56.28;
EngandWalesUK = 66.79/(56.28 + 3.15);

% set up and get data
%==========================================================================
spm_figure('GetWin','SI'); clf;
cd('C:\Users\karl\Dropbox\Coronavirus\Dashboard')


try
    
    % download data and write to CSV files
    %--------------------------------------------------------------------------
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newCasesBySpecimenDate&format=csv';
    writetable(webread(url,options),'cases.csv');
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newDeaths28DaysByDeathDate&format=csv';
    writetable(webread(url,options),'deaths.csv');
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=covidOccupiedMVBeds&format=csv';
    writetable(webread(url,options),'critical.csv');
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newPillarOneTwoTestsByPublishDate&format=csv';
    writetable(webread(url,options),'tests.csv');
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newAdmissions&format=csv';
    writetable(webread(url,options),'admissions.csv');
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=hospitalCases&format=csv';
    writetable(webread(url,options),'occupancy.csv');
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newOnsDeathsByRegistrationDate&format=csv';
    writetable(webread(url,options),'certified.csv');
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=uniqueCasePositivityBySpecimenDateRollingSum&format=csv';
    writetable(webread(url,options),'positivity.csv');
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newLFDTests&format=csv';
    writetable(webread(url,options),'lateralft.csv');
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=cumPeopleReceivingFirstDose&format=csv';
    writetable(webread(url,options),'vaccine.csv');
    
    url = 'https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/947572/COVID-19-transport-use-statistics.ods.ods';
    writetable(webread(url,options),'transport.csv');
    url = 'https://www.gstatic.com/covid19/mobility/2020_GB_Region_Mobility_Report.csv';
    tab = webread(url);
    writetable(tab(1:512,8:12),'mobility.csv');
    
    % url = 'https://www.ons.gov.uk/generator?uri=/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/articles/coronaviruscovid19infectionsinthecommunityinengland/december2020/b5e03a02&format=csv';
    % writetable(webread(url,options),'seropositive.csv');
    
    % https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/weekending11december2020
    url   = 'https://www.ons.gov.uk/generator?uri=/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/weekending11december2020/5b38d1d2&format=csv';
    writetable(webread(url,options),'place.csv');
    
    
    % https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-daily-deaths/
    url   = 'https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2020/12/COVID-19-total-announced-deaths-29-December-2020.xlsx';
    websave('ages.xlsx',url);
    
    disp('download successful')
    
catch
    
    disp('download failed')
    
end

% import data
%--------------------------------------------------------------------------
cases      = importdata('cases.csv');
deaths     = importdata('deaths.csv');
ccu        = importdata('critical.csv');
tests      = importdata('tests.csv');
certified  = importdata('certified.csv');
admissions = importdata('admissions.csv');
occupancy  = importdata('occupancy.csv');
positivity = importdata('positivity.csv');
lateralft  = importdata('lateralft.csv');
transport  = importdata('transport.csv');
mobility   = importdata('mobility.csv');
vaccine    = importdata('vaccine.csv');

serology   = importdata('seropositive.csv');
survey     = importdata('survey.csv');
symptoms   = importdata('symptoms.csv');
ratio      = importdata('ratio.csv');
place      = importdata('place.csv');

[num,txt]  = xlsread('ages.xlsx',5,'E16:KP23');

% created data structure
%--------------------------------------------------------------------------
Y(1).type = 'PCR cases (ONS)'; % daily PCR positive cases (by specimen)
Y(1).unit = 'number/day';
Y(1).U    = 2;
Y(1).date = datenum(cases.textdata(5:end,1),'yyyy-mm-dd');
Y(1).Y    = cases.data(4:end,1);
Y(1).h    = 0;
Y(1).lag  = 1;

Y(2).type = 'PCR tests (ONS)'; % daily PCR tests performed
Y(2).unit = 'number/day';
Y(2).U    = 6;
Y(2).date = datenum(tests.textdata(2:end,1),'yyyy-mm-dd');
Y(2).Y    = tests.data(:,1);
Y(2).h    = 0;
Y(2).lag  = 0;

Y(3).type = 'Virus/LFD tests (GOV)'; % newLFDTests (England)
Y(3).unit = 'number/day';
Y(3).U    = 24;
Y(3).date = datenum(lateralft.textdata(2:end,1),'yyyy-mm-dd');
Y(3).Y    = lateralft.data(:,1)*EnglandUK;
Y(3).h    = 0;
Y(3).lag  = 0;

Y(4).type = 'Prevalence (ONS)'; % number of people infected (England)
Y(4).unit = 'number';
Y(4).U    = 11;
Y(4).date = datenum(survey.textdata(2:end,1),'dd/mm/yyyy') - 2;
Y(4).Y    = survey.data(:,1)*EnglandUK;
Y(4).h    = 2;
Y(4).lag  = 0;

Y(5).type = 'Daily deaths (ONS: 28-days)'; % covid-related deaths (28 days)
Y(5).unit = 'number/day';
Y(5).U    = 1;
Y(5).date = datenum(deaths.textdata(8:end - 8,1),'yyyy-mm-dd');
Y(5).Y    = deaths.data(7:end - 8,1);
Y(5).h    = 0;
Y(5).lag  = 1;

Y(6).type = 'Certified deaths (ONS)'; % weekly covid related deaths
Y(6).unit = 'number';
Y(6).U    = 15;
Y(6).date = datenum(certified.textdata(2:end,1),'yyyy-mm-dd') - 10;
Y(6).Y    = certified.data(:,1)/7;
Y(6).h    = 0;
Y(6).lag  = 0;

Y(7).type = 'Admissions (ONS)'; % admissions to hospital
Y(7).unit = 'number';
Y(7).U    = 16;
Y(7).date = datenum(admissions.textdata(2:end,1),'yyyy-mm-dd');
Y(7).Y    = admissions.data(:,1);
Y(7).h    = 0;
Y(7).lag  = 0;

Y(8).type = 'Occupancy (ONS)'; % Hospital cases
Y(8).unit = 'number';
Y(8).U    = 27;
Y(8).date = datenum(occupancy.textdata(2:end,1),'yyyy-mm-dd');
Y(8).Y    = occupancy.data(:,1);
Y(8).h    = 0;
Y(8).lag  = 1;

Y(9).type = 'Ventilated patients (ONS)'; % CCU occupancy (mechanical)
Y(9).unit = 'number';
Y(9).U    = 3;
Y(9).date = datenum(ccu.textdata(2:end,1),'yyyy-mm-dd');
Y(9).Y    = ccu.data(:,1);
Y(9).h    = 2;
Y(9).lag  = 0;

Y(10).type = 'Seropositive (GOV)'; % percentage seropositive
Y(10).unit = 'percent';
Y(10).U    = 5;
Y(10).date = [datenum(serology.textdata(2:end,1),'dd/mm/yyyy') + 1; ...
    datenum(serology.textdata(2:end,1),'dd/mm/yyyy') + 2];
Y(10).Y    = [serology.data(:,2); serology.data(:,3)];
Y(10).h    = 0;
Y(10).lag  = 0;

Y(11).type = 'Symptoms (KCL)'; % number of people reporting symptoms (UK)
Y(11).unit = 'number';
Y(11).U    = 12;
Y(11).date = datenum(symptoms.textdata(2:end,1),'dd/mm/yyyy');
Y(11).Y    = symptoms.data(:,1);
Y(11).h    = 0;
Y(11).lag  = 1;

Y(12).type = 'R-ratio (WHO/GOV)'; % the production ratio
Y(12).unit = 'ratio';
Y(12).U    = 4;
Y(12).date = [datenum(ratio.textdata(2:end,1),'dd/mm/yyyy') - 7; ...
              datenum(ratio.textdata(2:end,1),'dd/mm/yyyy') - 8];
Y(12).Y    = [ratio.data(:,1); ratio.data(:,2)];
Y(12).h    = 0;
Y(12).lag  = 1;

Y(13).type = 'Transport (GOV)'; % cars (percent)
Y(13).unit = 'percent';
Y(13).U    = 13;
Y(13).Y    = transport.data(:,1)*100;
Y(13).date = datenum(transport.textdata(1 + (1:numel(Y(13).Y)),1),'dd-mm-yyyy');
Y(13).h    = 0;
Y(13).lag  = 0;

Y(14).type = 'Retail (Google)'; % retail and recreation (percent)
Y(14).unit = 'percent';
Y(14).U    = 14;
Y(14).date = datenum(mobility.textdata(2:end,1),'yyyy-mm-dd');
Y(14).Y    = mobility.data(:,1) + 100;
Y(14).h    = 0;
Y(14).lag  = 0;

Y(15).type = 'Hospital deaths (PHE)'; % hospital deaths
Y(15).unit = 'number';
Y(15).U    = 17;
Y(15).date = datenum(place.textdata(2:end,1),'dd mmm') - 1 - 365;
Y(15).Y    = place.data(:,4)*EngandWalesUK;
Y(15).h    = 0;
Y(15).lag  = 0;

Y(16).type = 'Hospital/Other deaths (PHE)'; % nonhospital deaths
Y(16).unit = 'number';
Y(16).U    = 18;
Y(16).date = datenum(place.textdata(2:end,1),'dd mmm') - 11 - 365;
Y(16).Y    = sum(place.data(:,1:3),2)*EngandWalesUK;
Y(16).h    = 0;
Y(16).lag  = 0;

Y(17).type = 'deaths > 60 (PHE)'; % deaths (English hospitals)
Y(17).unit = 'number';
Y(17).U    = 19;
Y(17).date = datenum(txt,'dd/mm/yyyy') - 1;
Y(17).Y    = sum(num(6:7,:))'*EnglandUK;
Y(17).h    = 0;
Y(17).lag  = 0;

Y(18).type = 'deaths </> 60 (PHE)'; % deaths (English hospitals)
Y(18).unit = 'number';
Y(18).U    = 20;
Y(18).date = datenum(txt,'dd/mm/yyyy') - 3;
Y(18).Y    = sum(num(3:5,:))'*EnglandUK;
Y(18).h    = 0;
Y(18).lag  = 0;

Y(19).type = 'PCR positivity (GOV)'; % positivity (England)
Y(19).unit = 'percent';
Y(19).U    = 23;
Y(19).date = datenum(positivity.textdata(2:end,1),'yyyy-mm-dd');
Y(19).Y    = positivity.data(:,1);
Y(19).h    = 0;
Y(19).lag  = 1;

Y(20).type = 'Vaccination (GOV)'; % New first dose
Y(20).unit = 'number';
Y(20).U    = 22;
Y(20).date = datenum(vaccine.textdata(2:end,1),'yyyy-mm-dd');
Y(20).Y    = vaccine.data(:,1);
Y(20).h    = 0;
Y(20).lag  = 0;

% normalise total deaths to ONS certified deaths
%--------------------------------------------------------------------------
d       = certified.data(:,1);
N       = sum(d(isfinite(d)));

d       = Y(15).Y;
n       = sum(d(isfinite(d)));
d       = Y(16).Y;
n       = sum(d(isfinite(d))) + n;
Y(15).Y = Y(15).Y*N/n;
Y(16).Y = Y(16).Y*N/n;

d       = Y(17).Y;
n       = sum(d(isfinite(d)));
d       = Y(18).Y;
n       = sum(d(isfinite(d))) + n;
Y(17).Y = Y(17).Y*N/n;
Y(18).Y = Y(18).Y*N/n;

% remove NANs, smooth and sort by date
%==========================================================================
M.date  = datestr(min(spm_vec(Y.date)),'dd-mm-yyyy');
Y(20)   = [];
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
pE.ho   = log([1 1]);           % coefficients for admissions
pC.ho   = [1 1]/8;              % prior variance
pE.hc   = log([1 1]);           % coefficients for hospital cases
pC.hc   = [1 1]/8;              % prior variance

pE.mo   = log([1.5,0.3]);       % coefficients for mobility
pC.mo   = [1 1]/8;              % prior variance
pE.wo   = log([1.5,0.3]);       % coefficients for workplace
pC.wo   = [1 1]/8;              % prior variance

pE.ag   = zeros(2,3);           % coefficients for age-related deaths
pC.ag   = ones(size(pE.ag));    % prior variance

% reporting lags
%--------------------------------------------------------------------------
lag([Y.U]) = [Y.lag];
pE.lag     = spm_zeros(lag);    % reporting delays
pC.lag     = lag;               % prior variance

% augment priors with fluctuations
%--------------------------------------------------------------------------
pE.vir     = zeros(1,8);        % transmission strength
pC.vir     = ones(1,8)/256;     % prior variance

pE.mob     = zeros(1,16);       % mobility
pC.mob     = ones(1,16)/8;      % prior variance

pE.pcr     = zeros(1,8);        % testing
pC.pcr     = ones(1,8)/64;      % prior variance

% shrink prior covariances for model comparison
%--------------------------------------------------------------------------
param = {'vir','mob','pcr'};
for f = 1:numel(param)
    if ~any(ismember(fluct,param{f}))
        pC.(param{f}) = pC.(param{f})*exp(-8);
    end
end


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
M.Nmax = 64;                   % maximum number of iterations
M.G    = @spm_SARS_gen;        % generative function
M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
M.pE   = pE;                   % prior expectations (parameters)
M.pC   = pC;                   % prior covariances  (parameters)
M.hE   = hE;                   % prior expectation  (log-precision)
M.hC   = 1/64;                 % prior covariances  (log-precision)
M.T    = Y;                    % data structure
U      = spm_vec(Y.U);         % outputs to model

% model inversion with Variational Laplace (Gauss Newton)
%==========================================================================
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,xY);

% save in DCM structure
%--------------------------------------------------------------------------
DCM.M  = M;
DCM.Ep = Ep;
DCM.Eh = Eh;
DCM.Cp = Cp;
DCM.F  = F;
DCM.Y  = YS;
DCM.U  = U;

% return if just DCM is required
%--------------------------------------------------------------------------
if nargout, return, end

% posterior predictions
%==========================================================================
spm_figure('GetWin','United Kingdom'); clf;
%--------------------------------------------------------------------------
M.T       = datenum('01-04-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
u         = [find(U == 1) find(U == 2) find(U == 3)];
[H,X,Z,R] = spm_SARS_gen(Ep,M,[1 2 3]);
spm_SARS_plot(H,X,YS(:,u),[1 2 3])

spm_figure('GetWin','outcomes (1)'); clf;
%--------------------------------------------------------------------------
j     = 0;
for i = 1:numel(Y)
    
    j = j + 1;
    subplot(4,2,j)
    spm_SARS_ci(Ep,Cp,YS(:,i),U(i),M);
    title(Y(i).type,'FontSize',14), ylabel(Y(i).unit)
    
    % add R = 1 and current date
    %----------------------------------------------------------------------
    if Y(i).U == 4
        plot(get(gca,'XLim'),[1,1],'-.r')
        plot(datenum(date)*[1,1],get(gca,'YLim'),'-.b')
        g  = get(gca,'YLim');
        plot(datenum(date)*[1,1],g,'-.b')
        set(gca,'YLim',[0 min(6,g(2))]), ylabel('ratio')
    end
    
    % hold plot
    %----------------------------------------------------------------------
    if ismember(Y(i).U,[6 17 19])
        j = j - 1; hold on
    end
    
    % new figure
    %----------------------------------------------------------------------
    if j == 8
        spm_figure('GetWin','outcomes (2)');
        j = 0;
    end
    
end

% infection fatality ratios (%)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes (3)');clf
j = 1;
subplot(4,2,j)
spm_SARS_ci(Ep,Cp,[],21,M);
ylabel('percent'),  title('Infection fatality ratio','FontSize',14)

% attack rate and herd immunity
%--------------------------------------------------------------------------
j = j + 1;
subplot(4,2,j)
spm_SARS_ci(Ep,Cp,[],25,M); hold on
spm_SARS_ci(Ep,Cp,[],26,M); hold off
ylabel('percent'),  title('Attack rate and immunity','FontSize',14)

% transmission risk
%--------------------------------------------------------------------------
j    = j + 1;
subplot(4,2,j)
plot([R.Ptrn]), spm_axis tight
title('Transmission risk','FontSize',14)
xlabel('days'),ylabel('probability')
hold on, plot([1,1]*size(DCM.Y,1),[0,1/2],':'), hold off, box off

% contact rate
%--------------------------------------------------------------------------
j    = j + 1;
subplot(4,2,j)
plot([R.Pout]), spm_axis tight
title('Contact rate','FontSize',14)
xlabel('days'),ylabel('probability')
hold on, plot([1,1]*size(DCM.Y,1),[0,1/2],':'), hold off, box off

j    = j + 1;
subplot(4,2,j)
plot(100 * [R.Pfat].*[R.Psev]), spm_axis tight
title('Symptom fatality ratio','FontSize',14)
xlabel('days'),ylabel('percent')
hold on, plot([1,1]*size(DCM.Y,1),[0,1/2],':'), hold off, box off

j    = j + 1;
subplot(4,2,j)
plot(100 * [R.Psen]), hold on, plot(100 * [R.Ptes]), spm_axis tight
title('Testing rate','FontSize',14)
xlabel('days'),ylabel('percent')
hold on, plot([1,1]*size(DCM.Y,1),[0,5],':'), hold off, box off
legend({'susceptible','infected'})

% save figures
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes (1)');
savefig(gcf,'Fig1')

spm_figure('GetWin','outcomes (2)');
savefig(gcf,'Fig2')

spm_figure('GetWin','outcomes (3)');
savefig(gcf,'Fig3')

spm_figure('GetWin','United Kingdom');
savefig(gcf,'Fig4')

% Table
%--------------------------------------------------------------------------
Tab = spm_COVID_table(Ep,Cp,M)

save('DCM_UK.mat','DCM')
cd('C:\Users\karl\Dropbox\Coronavirus')
save('DCM_UK.mat','DCM')

return


%% comparison of the effect of fluctuations
%==========================================================================

% Bayesian model reduction
%------------------------------------------------------------------------
fluct{1} = {'pcr'};
fluct{2} = {'vir','pcr'};
fluct{3} = {'mob','pcr'};
fluct{4} = {'vir','mob','pcr'};
fluct{5} = {};
fluct{6} = {'vir'};
fluct{7} = {'mob'};
fluct{8} = {'vir','mob'};

for i = 1:numel(fluct)
    DCM  = DEM_COVID_UK(fluct{i});
    F(i) = DCM.F
end

F = 1.0e+03 * [-7.4853   -7.1295   -7.1968   -6.9203   -7.5797   -7.1630   -7.2362   -6.9540]
F = F(:) - min(F);
P = spm_softmax(F(:));

spm_figure('GetWin','model comparison');clf
subplot(2,2,1), bar(F), axis square, box off
title('Model comparison','FontSize',14)
xlabel('model'), ylabel('log evidence')

subplot(2,2,2), bar(P), axis square, box off
title('Model comparison','FontSize',14)
xlabel('model'), ylabel('probability')




%% repeat with rapid loss of immunity
%==========================================================================

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
spm_SARS_plot(Z,X,DCM.YS,[2 1 3])
for j = 1:6
    subplot(3,2,j), hold on
    set(gca,'ColorOrderIndex',1);
end

M.CBT = CBT;
M.CBD = CBD;
[Z,X] = spm_SARS_gen(DCM.Ep,M,[2 1 3]);
spm_SARS_plot(Z,X,DCM.YS,[2 1 3])

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



%% fluctuations (adiabatic mean field approximation)
%==========================================================================
for f = 1:numel(fluct)
    
    % augment priors
    %----------------------------------------------------------------------
    pE.(fluct{f})   = zeros(1,16);           % add new prior expectation
    pC.(fluct{f})   =  ones(1,16);           % add new prior covariance
    
    % augment posteriors
    %----------------------------------------------------------------------
    i               = 1:size(Cp,1);          % number of parameters
    C               = Cp;                    % empirical prior covariance
    M.pE            = Ep;                    % empirical prior expectation
    M.pC            = spm_zeros(M.pC);       % fix expectations
    M.pE.(fluct{f}) = zeros(1,16);           % add new prior expectation
    M.pC.(fluct{f}) =  ones(1,16);           % add new prior covariance
    
    [Ep,Cp,Eh,Ff]   = spm_nlsi_GN(M,U,xY);   % new posterior expectation
    Cp(i,i)         = C;                     % new posterior covariance
    F               = F + Ff;                % free energy
    
    % save priors
    %----------------------------------------------------------------------
    M.pE   = pE;
    M.pC   = pC;
    
end





