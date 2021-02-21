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
% $Id: DEM_COVID_UK.m 8067 2021-02-21 16:15:48Z karl $

% DCM.F 06/02/2021: -1.8784e+04

% set up and preliminaries
%==========================================================================
% https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
% https://www.ndm.ox.ac.uk/covid-19/covid-19-infection-survey/results
% https://coronavirus.data.gov.uk/
% https://covid.joinzoe.com/data#levels-over-time
% https://www.gov.uk/guidance/the-r-number-in-the-uk#history
% https://www.gov.uk/government/statistics/transport-use-during-the-coronavirus-covid-19-pandemic
% https://www.google.com/covid19/mobility/

% mem = 256:  F = -1.0257e+04 
% mem = 2048: F = -1.0247e+04 

% web options
%--------------------------------------------------------------------------
options = weboptions('ContentType','table');
options.Timeout = 20;

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
    url = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=cumPeopleVaccinatedFirstDoseByVaccinationDate&format=csv';
    writetable(webread(url,options),'vaccine.csv');
    
    % get death by age
    %----------------------------------------------------------------------
    url   = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv';
    tab   = webread(url,options);
    age   = unique(tab(:,6));
    for r = 1:numel(age)
        j = find(ismember(tab(:,6),age(r,1)));
        ages(:,1)     = tab(j,1);
        ages(:,r + 1) = tab(j,7);
    end
    ages = renamevars(ages,(1:numel(age)) + 1,table2array(age));
    writetable(ages,'ages.csv')
    
    % get cases by age
    %----------------------------------------------------------------------
    url   = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newCasesBySpecimenDateAgeDemographics&format=csv';
    tab   = webread(url,options);
    age   = unique(tab(:,6));
    for r = 1:numel(age)
        j = find(ismember(tab(:,6),age(r,1)));
        agecases(:,1)     = tab(j,1);
        agecases(:,r + 1) = tab(j,7);
    end
    agecases = renamevars(agecases,(1:numel(age)) + 1,table2array(age));
    writetable(agecases,'agecases.csv')
    
    % mobility and transport
    %----------------------------------------------------------------------
    url = 'https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/947572/COVID-19-transport-use-statistics.ods.ods';
    writetable(webread(url,options),'transport.csv');
    url = 'https://www.gstatic.com/covid19/mobility/2020_GB_Region_Mobility_Report.csv';
    tab = webread(url);
    writetable(tab(1:512,9:12),'mobility.csv');
    

    disp('download successful')
    
catch
    
    disp('download failed')
    
end

% Files to be updated by hand
%--------------------------------------------------------------------------
if false
    url = 'https://www.ons.gov.uk/generator?uri=/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/weekending5february2021/ea5bdef1&format=csv';
    tab = webread(url,options);
    nd  = size(tab,1);
    dn  = (1:nd) + datenum(tab{end,1}) - nd;
    ds  = datestr(dn,'dd/mm/yyyy');
    for i = 1:nd
        tab(i,1) = {ds(i,:)};
    end
    writetable(tab,'place.csv');
    
    url = 'https://www.ons.gov.uk/generator?uri=/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/articles/coronaviruscovid19infectionsinthecommunityinengland/december2020/b5e03a02&format=csv';
    writetable(webread(url,options),'seropositive.csv');
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
ages       = importdata('ages.csv');
agecases   = importdata('agecases.csv');
place      = importdata('place.csv');

serology   = importdata('seropositive.csv');
survey     = importdata('survey.csv');
symptoms   = importdata('symptoms.csv');
ratio      = importdata('ratio.csv');


% create data structure
%--------------------------------------------------------------------------
Y(1).type = 'Positive virus tests (ONS)'; % daily PCR positive cases (by specimen)
Y(1).unit = 'number/day';
Y(1).U    = 2;
Y(1).date = datenum(cases.textdata(5:end,1),'yyyy-mm-dd');
Y(1).Y    = cases.data(4:end,1);
Y(1).h    = 0;
Y(1).lag  = 1;
Y(1).age  = 0;
Y(1).hold = 0;

Y(2).type = 'Virus tests (ONS)'; % daily PCR tests performed
Y(2).unit = 'number/day';
Y(2).U    = 6;
Y(2).date = datenum(tests.textdata(2:end,1),'yyyy-mm-dd');
Y(2).Y    = tests.data(:,1);
Y(2).h    = 0;
Y(2).lag  = 0;
Y(2).age  = 0;
Y(2).hold = 1;

Y(3).type = 'Virus/LFD tests (GOV)'; % newLFDTests (England)
Y(3).unit = 'number/day';
Y(3).U    = 24;
Y(3).date = datenum(lateralft.textdata(2:end,1),'yyyy-mm-dd');
Y(3).Y    = lateralft.data(:,1);
Y(3).h    = 0;
Y(3).lag  = 0;
Y(3).age  = 0;
Y(3).hold = 0;

Y(4).type = 'Prevalence (ONS)'; % number of people infected (England)
Y(4).unit = 'percent';
Y(4).U    = 11;
Y(4).date = datenum(survey.textdata(2:end,1),'dd/mm/yyyy');
Y(4).Y    = survey.data(:,1)*100;
Y(4).h    = 0;
Y(4).lag  = 1;
Y(4).age  = 0;
Y(4).hold = 0;

Y(5).type = 'Daily deaths (ONS: 28-days)'; % covid-related deaths (28 days)
Y(5).unit = 'number/day';
Y(5).U    = 1;
Y(5).date = datenum(deaths.textdata(5:end,1),'yyyy-mm-dd');
Y(5).Y    = deaths.data(4:end,1);
Y(5).h    = 2;
Y(5).lag  = 1;
Y(5).age  = 0;
Y(5).hold = 0;

Y(6).type = 'Certified deaths (ONS)'; % weekly covid related deaths
Y(6).unit = 'number';
Y(6).U    = 15;
Y(6).date = datenum(certified.textdata(2:end,1),'yyyy-mm-dd') - 10;
Y(6).Y    = certified.data(:,1)/7;
Y(6).h    = 0;
Y(6).lag  = 0;
Y(6).age  = 0;
Y(6).hold = 0;

Y(7).type = 'Admissions (ONS)'; % admissions to hospital
Y(7).unit = 'number';
Y(7).U    = 16;
Y(7).date = datenum(admissions.textdata(2:end,1),'yyyy-mm-dd');
Y(7).Y    = admissions.data(:,1);
Y(7).h    = 0;
Y(7).lag  = 0;
Y(7).age  = 0;
Y(7).hold = 0;

Y(8).type = 'Occupancy (ONS)'; % Hospital cases
Y(8).unit = 'number';
Y(8).U    = 27;
Y(8).date = datenum(occupancy.textdata(2:end,1),'yyyy-mm-dd');
Y(8).Y    = occupancy.data(:,1);
Y(8).h    = 0;
Y(8).lag  = 1;
Y(8).age  = 0;
Y(8).hold = 0;

Y(9).type = 'Ventilated patients (ONS)'; % CCU occupancy (mechanical)
Y(9).unit = 'number';
Y(9).U    = 3;
Y(9).date = datenum(ccu.textdata(2:end,1),'yyyy-mm-dd');
Y(9).Y    = ccu.data(:,1);
Y(9).h    = 2;
Y(9).lag  = 0;
Y(9).age  = 0;
Y(9).hold = 0;

Y(10).type = 'PCR positivity (GOV)'; % positivity (England)
Y(10).unit = 'percent';
Y(10).U    = 23;
Y(10).date = datenum(positivity.textdata(2:end,1),'yyyy-mm-dd');
Y(10).Y    = positivity.data(:,1);
Y(10).h    = 0;
Y(10).lag  = 1;
Y(10).age  = 0;
Y(10).hold = 0;

Y(11).type = 'Symptoms (KCL)'; % number of people reporting symptoms (UK)
Y(11).unit = 'number';
Y(11).U    = 12;
Y(11).date = datenum(symptoms.textdata(2:end,1),'dd/mm/yyyy');
Y(11).Y    = symptoms.data(:,1);
Y(11).h    = 0;
Y(11).lag  = 1;
Y(11).age  = 0;
Y(11).hold = 0;

Y(12).type = 'R-ratio (WHO/GOV)'; % the production ratio
Y(12).unit = 'ratio';
Y(12).U    = 4;
Y(12).date = [datenum(ratio.textdata(2:end,1),'dd/mm/yyyy') - 16; ...
              datenum(ratio.textdata(2:end,1),'dd/mm/yyyy') - 15];
Y(12).Y    = [ratio.data(:,1); ratio.data(:,2)];
Y(12).h    = 2;
Y(12).lag  = 1;
Y(12).age  = 0;
Y(12).hold = 0;

Y(13).type = 'Transport (GOV)'; % cars (percent)
Y(13).unit = 'percent';
Y(13).U    = 13;
Y(13).Y    = transport.data(:,1)*100;
Y(13).date = datenum(transport.textdata(1 + (1:numel(Y(13).Y)),1),'dd-mm-yyyy');
Y(13).h    = 0;
Y(13).lag  = 0;
Y(13).age  = 0;
Y(13).hold = 1;

Y(14).type = 'Mobility (GOV/Google)'; % retail and recreation (percent)
Y(14).unit = 'percent';
Y(14).U    = 14;
Y(14).date = datenum(mobility.textdata(2:end,1),'yyyy-mm-dd');
Y(14).Y    = mobility.data(:,1) + 100;
Y(14).h    = 0;
Y(14).lag  = 0;
Y(14).age  = 0;
Y(14).hold = 0;

% scaling for data from England
%--------------------------------------------------------------------------
EngWale    = sum(sum(place.data(1:end - 8,1:4),2));
UK         = sum(deaths.data(4:end,1));
EngWaleUK  = UK/EngWale;


Y(15).type = 'Hospital deaths (PHE)'; % hospital deaths
Y(15).unit = 'number';
Y(15).U    = 17;
Y(15).date = datenum(place.textdata(2:end - 8,1),'dd/mm/yyyy') - 1;
Y(15).Y    = place.data(1:end - 8,1)*EngWaleUK;
Y(15).h    = 0;
Y(15).lag  = 0;
Y(15).age  = 0;
Y(15).hold = 1;

Y(16).type = 'Hospital/Other deaths (PHE)'; % nonhospital deaths
Y(16).unit = 'number';
Y(16).U    = 18;
Y(16).date = datenum(place.textdata(2:end - 8,1),'dd/mm/yyyy') - 11;
Y(16).Y    = sum(place.data(1:end - 8,2:3),2)*EngWaleUK;
Y(16).h    = 0;
Y(16).lag  = 0;
Y(16).age  = 0;
Y(16).hold = 0;

Y(17).type = 'Seropositive (GOV)'; % percentage seropositive
Y(17).unit = 'percent';
Y(17).U    = 5;
Y(17).date = [datenum(serology.textdata(2:end,1),'dd/mm/yyyy') + 1; ...
              datenum(serology.textdata(2:end,1),'dd/mm/yyyy') + 2];
Y(17).Y    = [serology.data(:,2); serology.data(:,3)];
Y(17).h    = 2;
Y(17).lag  = 0;
Y(17).age  = 0;
Y(17).hold = 1;

Y(18).type = 'Ab+/Vaccination (GOV)'; % cumulative people with first dose
Y(18).unit = 'percent/millions';
Y(18).U    = 22;
Y(18).date = datenum(vaccine.textdata(2:end,1),'yyyy-mm-dd');
Y(18).Y    = vaccine.data(:,1)/1e6;
Y(18).h    = 0;
Y(18).lag  = 0;
Y(18).age  = 0;
Y(18).hold = 0;


% scaling for data from England
%--------------------------------------------------------------------------
England    = sum(sum(ages.data(:,[1 3 4 5 6 7:13 15 16:21]),2));
UK         = sum(deaths.data(4:end,1));
EnglandUK  = UK/England;


% age-specific data
%--------------------------------------------------------------------------
Y(19).type = 'deaths < 25 (PHE)'; % deaths (English hospitals)
Y(19).unit = 'number';
Y(19).U    = 1;
Y(19).date = datenum(ages.textdata(2:end,1),'yyyy-mm-dd');
Y(19).Y    = sum(ages.data(:,[1 3 4 5 6]),2)*EnglandUK;
Y(19).h    = 2;
Y(19).lag  = 0;
Y(19).age  = 1;
Y(19).hold = 1;

Y(20).type = 'deaths 25-65 (PHE)'; % deaths (English hospitals)
Y(20).unit = 'number';
Y(20).U    = 1;
Y(20).date = datenum(ages.textdata(2:end,1),'yyyy-mm-dd');
Y(20).Y    = sum(ages.data(:,[7:13 15]),2)*EnglandUK;
Y(20).h    = 2;
Y(20).lag  = 0;
Y(20).age  = 2;
Y(20).hold = 1;

Y(21).type = 'deaths 25-, 25-65, 65+ (PHE)'; % deaths (English hospitals)
Y(21).unit = 'number';
Y(21).U    = 1;
Y(21).date = datenum(ages.textdata(2:end,1),'yyyy-mm-dd');
Y(21).Y    = sum(ages.data(:,16:21),2)*EnglandUK;
Y(21).h    = 2;
Y(21).lag  = 0;
Y(21).age  = 3;
Y(21).hold = 0;


Y(22).type = 'cases < 25 (PHE)';   % cases (United Kingdom)
Y(22).unit = 'number';
Y(22).U    = 2;
Y(22).date = datenum(agecases.textdata(2:end,1),'yyyy-mm-dd');
Y(22).Y    = sum(agecases.data(:,[1 3 4 5 6]),2);
Y(22).h    = 0;
Y(22).lag  = 1;
Y(22).age  = 1;
Y(22).hold = 1;

Y(23).type = 'cases 25-65 (PHE)'; % cases (United Kingdom)
Y(23).unit = 'number';
Y(23).U    = 2;
Y(23).date = datenum(agecases.textdata(2:end,1),'yyyy-mm-dd');
Y(23).Y    = sum(agecases.data(:,[7:13 15]),2);
Y(23).h    = 0;
Y(23).lag  = 1;
Y(23).age  = 2;
Y(23).hold = 1;

Y(24).type = 'cases 25-, 25-65, 65+ (PHE)'; % cases (United Kingdom)
Y(24).unit = 'number';
Y(24).U    = 2;
Y(24).date = datenum(agecases.textdata(2:end,1),'yyyy-mm-dd');
Y(24).Y    = sum(agecases.data(:,16:21),2);
Y(24).h    = 0;
Y(24).lag  = 1;
Y(24).age  = 3;
Y(24).hold = 0;

% remove location dependent data
%--------------------------------------------------------------------------
% Y([15 16]) = [];


% population sizes
%--------------------------------------------------------------------------
% {'00_04','05_09','10_14','15_19','20_24','25_29','30_34','35_39',...
%  '40_44','45_49','50_54','55_59','60_64','65_69','70_74','75_79', ...
%  '80_84','85_89','90+'};
N  = [3.86, 4.15, 3.95, 3.66, 4.15, 4.51, 4.50, 4.40, 4.02, 4.4, 4.66, ...
      4.41, 3.76, 3.37, 3.32, 2.33, 1.72, 1.04, 0.61];
N  = [sum(N(1:5)) sum(N(6:13)) sum(N(14:end))]; % age bands
nN = numel(N);

% remove NANs, smooth and sort by date
%==========================================================================
M.date  = datestr(min(spm_vec(Y.date)),'dd-mm-yyyy');
M.date  = '01-02-2020';
[Y,YS]  = spm_COVID_Y(Y,M.date);

% data structure with vectorised data and covariance components
%--------------------------------------------------------------------------
xY.y    = spm_vec(Y.Y);
xY.Q    = spm_Ce([Y.n]);
hE      = spm_vec(Y.h);

% get and set priors
%--------------------------------------------------------------------------
[pE,pC] = spm_SARS_priors(nN);
pE.N    = log(N(:));
pC.N    = spm_zeros(pE.N);

% coefficients for likelihood model
%--------------------------------------------------------------------------
% pE.sy = log([1 1;1 1; 1 1]); % coefficients for symptoms
% pC.sy = [1 1;1 1; 1 1]/8;    % prior variance

pE.dc   = log([1 1/16]);       % coefficients for death (28 days)
pC.dc   = [1 1]/8;             % prior variance
pE.ho   = log([8 1/8]);        % coefficients for admissions
pC.ho   = [1 1]/8;             % prior variance
pE.hc   = log([8 1/4]);        % coefficients for hospital cases
pC.hc   = [1 1]/8;             % prior variance
pE.mo   = log([1 1]);          % coefficients for mobility
pC.mo   = [1 1]/8;             % prior variance
pE.wo   = log([1 1]);          % coefficients for workplace
pC.wo   = [1 1]/8;             % prior variance

pE.ps   = log(1);              % coefficient for positivity estimate
pC.ps   = 1/8;                 % prior variance

% reporting lags
%--------------------------------------------------------------------------
lag([Y.U]) = [Y.lag];

pE.lag  = spm_zeros(lag);      % reporting delays
pC.lag  = lag;                 % prior variance

% augment priors with fluctuations
%--------------------------------------------------------------------------
pE.tra  = zeros(1,8);          % transmission strength
pC.tra  = ones(1,8)/256;       % prior variance

pE.pcr  = zeros(1,8);          % testing
pC.pcr  = ones(1,8)/8;       % prior variance

pE.mob  = zeros(1,16);         % mobility
pC.mob  = ones(1,16)/8;      % prior variance

% model specification
%==========================================================================
M.Nmax = 32;                   % maximum number of iterations
M.G    = @spm_SARS_gen;        % generative function
M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
M.pE   = pE;                   % prior expectations (parameters)
M.pC   = pC;                   % prior covariances  (parameters)
M.hE   = hE;                   % prior expectation  (log-precision)
M.hC   = 1/512;                % prior covariances  (log-precision)
M.T    = Y;                    % data structure

U      = [Y.U];                % outputs to model
A      = [Y.age];              % age bands

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
DCM.A  = A;

% return if just DCM is required
%--------------------------------------------------------------------------
if nargout, return, end

% posterior predictions
%==========================================================================
spm_figure('GetWin','United Kingdom'); clf;
%--------------------------------------------------------------------------
M.T       = datenum('01-04-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
u         = [find(U == 1,1) find(U == 2,1) find(U == 3,1)];
[H,X,~,R] = spm_SARS_gen(Ep,M,[1 2 3]);
spm_SARS_plot(H,X,YS(:,u),[1 2 3])

spm_figure('GetWin','outcomes (1)');
%--------------------------------------------------------------------------
j     = 0;
for i = 1:numel(Y)
    
    j = j + 1;
    subplot(4,2,j)
    spm_SARS_ci(Ep,Cp,YS(:,i),U(i),M,[],A(i));
    title(Y(i).type,'FontSize',14), ylabel(Y(i).unit)
    
    % add R = 1 and current date
    %----------------------------------------------------------------------
    if Y(i).U == 4
        plot(get(gca,'XLim'),[1,1],'-.r')
        plot(datenum(date)*[1,1],get(gca,'YLim'),'-.b')
        set(gca,'YLim',[0 5]), ylabel('ratio')
    end
    
    % hold plot
    %----------------------------------------------------------------------
    if Y(i).hold
        j = j - 1; hold on
    end
    
    % new figure
    %----------------------------------------------------------------------
    if j == 8
        spm_figure('GetWin','outcomes (2)');
        j = 0;
    end
    
end

% time varying parameters
%==========================================================================

% infection fatality ratios (%)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes (3)');
j = 1;
subplot(4,2,j)
for i = 1:numel(N)
    spm_SARS_ci(Ep,Cp,[],21,M,[],i); hold on
end
ylabel('percent'), title('Infection fatality ratio','FontSize',14)

% transmission risk
%--------------------------------------------------------------------------
j    = j + 1;
subplot(4,2,j), hold on
plot([R{1}.Ptrn]), spm_axis tight
title('Transmission risk','FontSize',14)
xlabel('days'),ylabel('probability')
hold on, plot([1,1]*size(DCM.Y,1),[0,1/2],':'), box off

% contact rate
%--------------------------------------------------------------------------
j    = j + 1;
subplot(4,2,j), hold on
for i = 1:numel(R)
    plot([R{i}.Pout])
end
spm_axis tight
title('Contact rate','FontSize',14)
xlabel('days'),ylabel('probability')
hold on, plot([1,1]*size(DCM.Y,1),[0,1/2],':'), box off


% case fatality ratio
%--------------------------------------------------------------------------
j    = j + 1;
subplot(4,2,j), hold on
for i = 1:numel(R)
    plot(100 * [R{i}.Pfat])
end
spm_axis tight
title('Fatality risk | ARDS','FontSize',14)
xlabel('days'),ylabel('percent')
hold on, plot([1,1]*size(DCM.Y,1),[0,1/2],':'), box off
legend({'< 24yrs','25-64yrs','> 64yrs'})

% testing rates
%--------------------------------------------------------------------------
% j    = j + 1;
% subplot(4,2,j), hold on
% for i = 1:numel(R)
%     plot(100 * [R{i}.Psen]), plot(100 * [R{i}.Ptes])
% end
% spm_axis tight
% title('Testing rates (un/infected)','FontSize',14)
% xlabel('days'),ylabel('percent')
% plot([1,1]*size(DCM.Y,1),[0,5],':'), box off


%% long-term forecasts
%==========================================================================
spm_figure('GetWin','outcomes (3)');

Ep  = DCM.Ep;
Cp  = DCM.Cp;
M   = DCM.M;
M.T = datenum('01-11-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
t   = (1:M.T) + datenum(M.date,'dd-mm-yyyy');

% infection fatality ratios (%)
%--------------------------------------------------------------------------
subplot(2,1,2)
spm_SARS_ci(Ep,Cp,[],19,M); hold on
spm_SARS_ci(Ep,Cp,[],20,M); hold on
ylabel('cases per 100,000'), title('Incidence per 100,000','FontSize',14)
plot(datenum(date)*[1,1],get(gca,'YLim'),':')
legend({'CI per day','actual cases per day','CI per week','confirmed cases per week'})

% switch windows
%--------------------------------------------------------------------------
spm_figure('GetWin','long-term (1)');

% fatalities
%--------------------------------------------------------------------------
subplot(2,1,1)

i   = find(DCM.U == 1,1);  D = DCM.Y(:,i); spm_SARS_ci(Ep,Cp,D,1,M);  hold on
i   = find(DCM.U == 15,1); D = DCM.Y(:,i); spm_SARS_ci(Ep,Cp,D,15,M); hold on

plot(datenum(date,'dd-mm-yyyy')*[1,1],get(gca,'YLim'),':k')
ylabel('number per day'), title('Daily deaths','FontSize',14)
legend({'CI 28-day','PCR test within 28 days','ONS','CI certified','certified deaths'})
legend boxoff
drawnow

% lockdown and mobility
%--------------------------------------------------------------------------
subplot(2,1,2)
i       = find(DCM.U == 14,1); D = DCM.Y(:,i);
q       = spm_SARS_gen(Ep,M,14);

[~,~,q] = spm_SARS_ci(Ep,Cp,D,14,M); hold on

% thresholds
%--------------------------------------------------------------------------
u1   = datenum('10-May-2020','dd-mmm-yyyy') - t(1) + 1;
u2   = datenum('20-Aug-2020','dd-mmm-yyyy') - t(1) + 1;
u3   = datenum('10-Mar-2021','dd-mmm-yyyy') - t(1) + 1;
U    = sort([0 q(u1) q(u2) q(u3)]);
dstr = datestr(t,'dd-mmm');

% loop over levels
%==========================================================================
for i = 1:numel(U)
    
    % intervals for this level
    %----------------------------------------------------------------------
    if i == 1
        j  = find(q <= U(i + 1));
    elseif i == numel(U)
        j  = find(q >= U(i));
    else
        j  = find(q >= U(i) & q <= U(i + 1));
    end
    
    % Timeline
    %----------------------------------------------------------------------
    for k = 1:numel(j)
        try
            fill(t(j(k) + [0 1 1 0]),[0 0 1 1]*32,'r', ...
                'FaceAlpha',(numel(U) - i)/16,'Edgecolor','none')
        end
    end
    
    % label level crossings
    %----------------------------------------------------------------------
    if i <numel(U)
        j = find((q(1:end - 1) <= U(i + 1)) & (q(2:end) > U(i + 1)));
    else
        j = [];
    end
    for k = 1:numel(j)
        text(t(j(k)),i*8,dstr(j(k),:),'Color','k','FontSize',10)
    end
    
    % plot levels
    %----------------------------------------------------------------------
    plot([t(1),t(end)],U(i)*[1,1],':r')
    
end

ylabel('percent'),  title('Mobility and lockdown','FontSize',14)
legend({'credible interval','mobility (%)'}), legend boxoff
drawnow


% prevalence and reproduction ratio
%--------------------------------------------------------------------------
spm_figure('GetWin','long-term (2)');
subplot(2,1,1)
M.T = datenum('01-11-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
t   = (1:M.T) + datenum(M.date,'dd-mm-yyyy');

i   = find(DCM.U == 4,1);
Rt  = DCM.Y(:,i);
            spm_SARS_ci(Ep,Cp,[],11,M); hold on
[~,~,q,c] = spm_SARS_ci(Ep,Cp,Rt,4 ,M); hold on

j   = find(t == datenum(date));
q   = q(j);
d   = sqrt(c{1}(j,j))*1.64;
str = sprintf('Prevalence and reproduction ratio (%s): R = %.2f (CI %.2f to %.2f)',datestr(date,'dd-mmm-yy'),q,q - d,q + d);


% add R = 1 and current dateline
%--------------------------------------------------------------------------
plot(get(gca,'XLim'),[1,1],':k')
plot(datenum(date,'dd-mm-yyyy')*[1,1],get(gca,'YLim'),':k')
set(gca,'YLim',[0 5]), ylabel('ratio / percent')
ylabel('percent / millions'),  title(str,'FontSize',14)

legend({'CI prevalence','Prevalence (%)','CI R-number','R DCM','R SPI-M'})
legend boxoff
drawnow

% attack rate, herd immunity and herd immunity threshold
%--------------------------------------------------------------------------
subplot(2,1,2)
spm_SARS_ci(Ep,Cp,[],25,M); hold on
spm_SARS_ci(Ep,Cp,[],26,M); hold on

[H,~,~,R] = spm_SARS_gen(Ep,M,[4 22 26]);
i         = 1:32;                           % pre-pandemic period
TRN       = [R{1}.Ptrn];                    % transmission risk
R0        = mean(H(i,1));                   % basic reproduction ratio    
RT        = R0*TRN(:)/mean(TRN(i));         % effective reproduction ratio
HIT       = 100 * (1 - 1./RT);              % herd immunity threshold
VAC       = H(:,2);                         % number of people vaccinated
i         = find(H(:,3) > HIT,1);           % date threshold reached
i         = min([i,M.T]);

q   = Ep.vac;
d   = spm_unvec(diag(Cp),Ep);
d   = sqrt(d.vac)*1.64;
qE  = 100*exp(q);
qL  = 100*exp(q - d);
qU  = 100*exp(q + d);
str = sprintf('Attack rate and immunity: vaccine efficacy %.1f%s (CI %.1f to %.1f)',qE,'%',qL,qU);
leg = sprintf('%s (HIT: %.1f%s)',datestr(t(i),'dd-mmm-yy'),HIT(i),'%');

plot(t,HIT,t,VAC), hold on
plot(t(i)*[1,1],[0,100],':k'), set(gca,'YLim',[0,100])
ylabel('percent / millions'),  title(str,'FontSize',14)
legend({'CI','Attack rate','CI','Herd immunity','Herd immunity threshold','Numer of vaccinations'})
legend boxoff
text(t(i),8,leg,'FontSize',10), drawnow


%% save figures
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes (1)');
savefig(gcf,'Fig1')

spm_figure('GetWin','outcomes (2)');
savefig(gcf,'Fig2')

spm_figure('GetWin','outcomes (3)');
savefig(gcf,'Fig3')

spm_figure('GetWin','United Kingdom');
savefig(gcf,'Fig4')

spm_figure('GetWin','long-term (2)');
savefig(gcf,'Fig5')

spm_figure('GetWin','long-term (1)');
savefig(gcf,'Fig6')

% Table
%--------------------------------------------------------------------------
Tab = spm_COVID_table(Ep,Cp,M)

save('DCM_UK.mat','DCM')
cd('C:\Users\karl\Dropbox\Coronavirus')
save('DCM_UK.mat','DCM')

return



%% NOTES

% fluctuations (adiabatic mean field approximation)
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





