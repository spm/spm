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
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% set up and preliminaries
%==========================================================================
% https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
% https://www.ndm.ox.ac.uk/covid-19/covid-19-infection-survey/results
% https://coronavirus.data.gov.uk/
% https://covid.joinzoe.com/data#levels-over-time
% https://www.gov.uk/guidance/the-r-number-in-the-uk#history
% https://www.gov.uk/government/statistics/transport-use-during-the-coronavirus-covid-19-pandemic
% https://www.google.com/covid19/mobility/


% set up and get data
%==========================================================================
spm_figure('GetWin','SI'); clf;
cd('C:\Users\karl\Dropbox\Coronavirus\Dashboard')

% Files to be updated by hand
%--------------------------------------------------------------------------
% url = 'https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021/20210414covid19infectionsurveydatasets.xlsx'
% tab = webread(url);
% url = 'https://www.ons.gov.uk/generator?uri=/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/weekending19february2021/8714ef2a&format=csv';
% writetable(webread(url,options),'place.csv');


%% web options
%--------------------------------------------------------------------------
options = weboptions('ContentType','table');
options.Timeout = 20;

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

% get death by age (England)
%--------------------------------------------------------------------------
url   = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv';
tab   = webread(url,options);
vnames = tab.Properties.VariableNames;
aa     = find(ismember(vnames,'age'));
ad     = find(ismember(vnames,'date'));
an     = find(ismember(vnames,'deaths'));
age   = unique(tab(:,aa));
for r = 1:numel(age)
    j = find(ismember(tab(:,aa),age(r,1)));
    agedeaths(:,1)     = tab(j,ad);
    agedeaths(:,r + 1) = tab(j,an);
end
agedeaths = renamevars(agedeaths,(1:numel(age)) + 1,table2array(age));
writetable(agedeaths,'agedeaths.csv')

% get cases by age (UK)
%----------------------------------------------------------------------
url   = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDateAgeDemographics&format=csv';
tab   = webread(url,options);
vnames = tab.Properties.VariableNames;
aa     = find(ismember(vnames,'age'));
ad     = find(ismember(vnames,'date'));
an     = find(ismember(vnames,'cases'));
age   = unique(tab(:,aa));
for r = 1:numel(age)
    j = find(ismember(tab(:,aa),age(r,1)));
    agecases(:,1)     = tab(j,ad);
    agecases(:,r + 1) = tab(j,an);
end
agecases = renamevars(agecases,(1:numel(age)) + 1,table2array(age));
writetable(agecases,'agecases.csv')

% mobility and transport
%--------------------------------------------------------------------------
url = 'https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/947572/COVID-19-transport-use-statistics.ods.ods';
writetable(webread(url,options),'transport.csv');

ndy = 321;
url = 'https://www.gstatic.com/covid19/mobility/2020_GB_Region_Mobility_Report.csv';
tab = webread(url);
writetable(tab(1:ndy,9:12),'mobility20.csv');

ndy = datenum(date) - datenum(datestr('01/01/2021','dd/mm/yyyy'));
url = 'https://www.gstatic.com/covid19/mobility/2021_GB_Region_Mobility_Report.csv';
tab = webread(url);
writetable(tab(1:ndy,9:12),'mobility21.csv');

% population sizes (millions)
%--------------------------------------------------------------------------
% {'00_04','05_09','10_14','15_19','20_24','25_29','30_34','35_39',...
%  '40_44','45_49','50_54','55_59','60_64','65_69','70_74','75_79', ...
%  '80_84','85_89','90+'};
N  = [3.86, 4.15, 3.95, 3.66, 4.15, 4.51, 4.50, 4.40, 4.02, 4.4, 4.66, ...
      4.41, 3.76, 3.37, 3.32, 2.33, 1.72, 1.04, 0.61];

% ONS age bands
%--------------------------------------------------------------------------
ons{1} = [sum(N(4:5))];
ons{2} = [sum(N(6:7)) sum(N(8:10)) sum(N(11:12)) N(13)];
ons{3} = [N(14:16) sum(N(17:end))];
vac    = ons{1}'/sum(N(1:5));
ons{1} = ons{1}'/sum(ons{1});
ons{2} = ons{2}'/sum(ons{2});
ons{3} = ons{3}'/sum(ons{3});

% DCM age bands
%--------------------------------------------------------------------------
N   = [sum(N(1:5)) sum(N(6:13)) sum(N(14:end))];
nN  = numel(N);




%% import data
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
mobility20 = importdata('mobility20.csv');
mobility21 = importdata('mobility21.csv');
vaccine    = importdata('vaccineage.csv');
agedeaths  = importdata('agedeaths.csv');
agecases   = importdata('agecases.csv');

serology   = importdata('serology.csv');
place      = importdata('place.csv');
survey     = importdata('survey.csv');
surveyage  = importdata('surveyage.csv');
symptoms   = importdata('symptoms.csv');
ratio      = importdata('ratio.csv');

d          = find(ismember(cases.textdata(1,1:end),'date'));

% create data structure
%--------------------------------------------------------------------------
Y(1).type = 'Positive virus tests (ONS)'; % daily PCR positive cases (by specimen)
Y(1).unit = 'number/day';
Y(1).U    = 2;
Y(1).date = datenum(cases.textdata(5:end,d),'yyyy-mm-dd');
Y(1).Y    = cases.data(4:end,1);
Y(1).h    = 0;
Y(1).lag  = 1;
Y(1).age  = 0;
Y(1).hold = 0;

Y(2).type = 'Virus tests (ONS)'; % daily PCR tests performed
Y(2).unit = 'number/day';
Y(2).U    = 6;
Y(2).date = datenum(tests.textdata(2:end,d),'yyyy-mm-dd');
Y(2).Y    = tests.data(:,1);
Y(2).h    = 0;
Y(2).lag  = 0;
Y(2).age  = 0;
Y(2).hold = 1;

Y(3).type = 'Virus/LFD tests (GOV)'; % newLFDTests (England)
Y(3).unit = 'number/day';
Y(3).U    = 24;
Y(3).date = datenum(lateralft.textdata(2:end,d),'yyyy-mm-dd');
Y(3).Y    = lateralft.data(:,1);
Y(3).h    = 0;
Y(3).lag  = 0;
Y(3).age  = 0;
Y(3).hold = 0;

Y(4).type = 'Prevalence (ONS)'; % number of people infected (England)
Y(4).unit = 'percent';
Y(4).U    = 11;
Y(4).date = datenum(survey.textdata(2:end,1),'dd/mm/yyyy') - 7;
Y(4).Y    = survey.data(:,1)*100;
Y(4).h    = 0;
Y(4).lag  = 1;
Y(4).age  = 0;
Y(4).hold = 0;

Y(5).type = 'Daily deaths (ONS: 28-days)'; % covid-related deaths (28 days)
Y(5).unit = 'number/day';
Y(5).U    = 1;
Y(5).date = datenum(deaths.textdata(5:end,d),'yyyy-mm-dd');
Y(5).Y    = deaths.data(4:end,1);
Y(5).h    = 2;
Y(5).lag  = 1;
Y(5).age  = 0;
Y(5).hold = 0;

Y(6).type = 'Certified deaths (ONS)'; % weekly covid related deaths
Y(6).unit = 'number';
Y(6).U    = 15;
Y(6).date = datenum(certified.textdata(2:end,d),'yyyy-mm-dd') - 10;
Y(6).Y    = certified.data(:,1)/7;
Y(6).h    = 2;
Y(6).lag  = 0;
Y(6).age  = 0;
Y(6).hold = 0;

Y(7).type = 'Admissions (ONS)'; % admissions to hospital
Y(7).unit = 'number';
Y(7).U    = 16;
Y(7).date = datenum(admissions.textdata(2:end,d),'yyyy-mm-dd');
Y(7).Y    = admissions.data(:,1);
Y(7).h    = 0;
Y(7).lag  = 0;
Y(7).age  = 0;
Y(7).hold = 0;

Y(8).type = 'Occupancy (ONS)'; % Hospital cases
Y(8).unit = 'number';
Y(8).U    = 27;
Y(8).date = datenum(occupancy.textdata(2:end,d),'yyyy-mm-dd');
Y(8).Y    = occupancy.data(:,1);
Y(8).h    = 0;
Y(8).lag  = 1;
Y(8).age  = 0;
Y(8).hold = 0;

Y(9).type = 'Ventilated patients (ONS)'; % CCU occupancy (mechanical)
Y(9).unit = 'number';
Y(9).U    = 3;
Y(9).date = datenum(ccu.textdata(2:end,d),'yyyy-mm-dd');
Y(9).Y    = ccu.data(:,1);
Y(9).h    = 2;
Y(9).lag  = 0;
Y(9).age  = 0;
Y(9).hold = 0;

Y(10).type = 'PCR positivity (GOV)'; % positivity (England)
Y(10).unit = 'percent';
Y(10).U    = 23;
Y(10).date = datenum(positivity.textdata(2:end,d),'yyyy-mm-dd');
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
Y(14).date = [datenum(mobility20.textdata(2:end,1),'yyyy-mm-dd') ;
    datenum(mobility21.textdata(2:end,1),'yyyy-mm-dd')];
Y(14).Y    = [mobility20.data(:,1); mobility21.data(:,1)] + 100;
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
Y(15).date = datenum(place.textdata(2:end - 8,1),'dd/mm/yyyy');
Y(15).Y    = place.data(1:end - 8,4)*EngWaleUK;
Y(15).h    = 0;
Y(15).lag  = 1;
Y(15).age  = 0;
Y(15).hold = 1;

Y(16).type = 'Hospital/Other deaths (PHE)'; % nonhospital deaths
Y(16).unit = 'number';
Y(16).U    = 18;
Y(16).date = datenum(place.textdata(2:end - 8,1),'dd/mm/yyyy');
Y(16).Y    = sum(place.data(1:end - 8,1:3),2)*EngWaleUK;
Y(16).h    = 0;
Y(16).lag  = 1;
Y(16).age  = 0;
Y(16).hold = 0;

% scaling for data from England
%--------------------------------------------------------------------------
ig1        = [1 3 4 5 6];                                % <25
ig2        = [7:13 15];                                  % 25-65
ig3        = [16:21];                                    % >65
ig         = [ig1 ig2 ig3];
England    = sum(sum(agedeaths.data(:,ig),2));
UK         = sum(deaths.data(4:end,1));
EnglandUK  = UK/England;

% age-specific data
%--------------------------------------------------------------------------
j          = find(~ismember(serology.textdata(1,2:end),''));
Y(17).type = 'Seropositive < 25 (PHE)'; % percent antibody positive (England)
Y(17).unit = 'percent';
Y(17).U    = 5;
Y(17).date = datenum(serology.textdata(2:end,1),'dd/mm/yyyy');
Y(17).Y    = serology.data(:,j(1))*ons{1}*100;
Y(17).h    = 2;
Y(17).lag  = 0;
Y(17).age  = 1;
Y(17).hold = 1;

Y(18).type = 'Seropositive 25-65 (PHE)'; % percent antibody positive (England)
Y(18).unit = 'percent';
Y(18).U    = 5;
Y(18).date = datenum(serology.textdata(2:end,1),'dd/mm/yyyy');
Y(18).Y    = serology.data(:,j(2:5))*ons{2}*100;
Y(18).h    = 2;
Y(18).lag  = 0;
Y(18).age  = 2;
Y(18).hold = 1;

Y(19).type = 'Seropositive 25-, 25-65, 65+ (PHE)'; % percent antibody positive (England)
Y(19).unit = 'percent';
Y(19).U    = 5;
Y(19).date = datenum(serology.textdata(2:end,1),'dd/mm/yyyy');
Y(19).Y    = serology.data(:,j(6:9))*ons{3}*100;
Y(19).h    = 2;
Y(19).lag  = 0;
Y(19).age  = 3;
Y(19).hold = 0;

j          = find(~ismember(vaccine.textdata(1,2:end),''));
Y(20).type = 'First dose < 25 (PHE)'; % percent vaccinated (England)
Y(20).unit = 'percent';
Y(20).U    = 22;
Y(20).date = datenum(vaccine.textdata(2:end,1),'dd/mm/yyyy');
Y(20).Y    = vaccine.data(:,j(1))*vac*100;
Y(20).h    = 2;
Y(20).lag  = 0;
Y(20).age  = 1;
Y(20).hold = 1;

Y(21).type = 'First dose 25-65 (PHE)'; % percent vaccinated (England)
Y(21).unit = 'percent';
Y(21).U    = 22;
Y(21).date = datenum(vaccine.textdata(2:end,1),'dd/mm/yyyy');
Y(21).Y    = vaccine.data(:,j(2:5))*ons{2}*100;
Y(21).h    = 2;
Y(21).lag  = 0;
Y(21).age  = 2;
Y(21).hold = 1;

Y(22).type = 'First dose 25-, 25-65, 65+ (PHE)'; % percent vaccinated (England)
Y(22).unit = 'percent';
Y(22).U    = 22;
Y(22).date = datenum(vaccine.textdata(2:end,1),'dd/mm/yyyy');
Y(22).Y    = vaccine.data(:,j(6:9))*ons{3}*100;
Y(22).h    = 2;
Y(22).lag  = 0;
Y(22).age  = 3;
Y(22).hold = 0;


Y(23).type = 'Deaths < 25 (PHE)'; % deaths (English hospitals)
Y(23).unit = 'number';
Y(23).U    = 1;
Y(23).date = datenum(agedeaths.textdata(2:end,1),'yyyy-mm-dd');
Y(23).Y    = sum(agedeaths.data(:,ig1),2)*EnglandUK;
Y(23).h    = 2;
Y(23).lag  = 0;
Y(23).age  = 1;
Y(23).hold = 1;

Y(24).type = 'Deaths 25-65 (PHE)'; % deaths (English hospitals)
Y(24).unit = 'number';
Y(24).U    = 1;
Y(24).date = datenum(agedeaths.textdata(2:end,1),'yyyy-mm-dd');
Y(24).Y    = sum(agedeaths.data(:,ig2),2)*EnglandUK;
Y(24).h    = 2;
Y(24).lag  = 0;
Y(24).age  = 2;
Y(24).hold = 1;

Y(25).type = 'Deaths 25-, 25-65, 65+ (PHE)'; % deaths (English hospitals)
Y(25).unit = 'number';
Y(25).U    = 1;
Y(25).date = datenum(agedeaths.textdata(2:end,1),'yyyy-mm-dd');
Y(25).Y    = sum(agedeaths.data(:,ig3),2)*EnglandUK;
Y(25).h    = 2;
Y(25).lag  = 0;
Y(25).age  = 3;
Y(25).hold = 0;


Y(26).type = 'PCR cases < 25 (PHE)';  % PCR notifications (United Kingdom)
Y(26).unit = 'number';
Y(26).U    = 2;
Y(26).date = datenum(agecases.textdata(2:end,1),'yyyy-mm-dd');
Y(26).Y    = sum(agecases.data(:,ig1),2);
Y(26).h    = 0;
Y(26).lag  = 1;
Y(26).age  = 1;
Y(26).hold = 1;

Y(27).type = 'PCR cases 25-65 (PHE)'; % PCR notifications  (United Kingdom)
Y(27).unit = 'number';
Y(27).U    = 2;
Y(27).date = datenum(agecases.textdata(2:end,1),'yyyy-mm-dd');
Y(27).Y    = sum(agecases.data(:,ig2),2);
Y(27).h    = 0;
Y(27).lag  = 1;
Y(27).age  = 2;
Y(27).hold = 1;

Y(28).type = 'PCR cases 25-, 25-65, 65+ (PHE)'; % PCR notifications  (United Kingdom)
Y(28).unit = 'number';
Y(28).U    = 2;
Y(28).date = datenum(agecases.textdata(2:end,1),'yyyy-mm-dd');
Y(28).Y    = sum(agecases.data(:,ig3),2);
Y(28).h    = 0;
Y(28).lag  = 1;
Y(28).age  = 3;
Y(28).hold = 0;

j          = find(~ismember(surveyage.textdata(1,2:end),''));
cj         = [mean(surveyage.data(end,j(1:3)),2);
              mean(surveyage.data(end,j(4:5)),2);
              mean(surveyage.data(end,j(6:7)),2)];
cj         = survey.data(end,1)/( (N*cj)/sum(N) );
Y(29).type = 'Prevalence < 25 (PHE)';  % Estimated positivity (England)
Y(29).unit = 'percent';
Y(29).U    = 11;
Y(29).date = datenum(surveyage.textdata(2:end,1),'dd/mm/yyyy');
Y(29).Y    = mean(surveyage.data(:,j(1:3)),2)*100*cj;
Y(29).h    = 0;
Y(29).lag  = 0;
Y(29).age  = 1;
Y(29).hold = 1;

Y(30).type = 'Prevalence 25-65 (PHE)'; % Estimated positivity (England)
Y(30).unit = 'percent';
Y(30).U    = 11;
Y(30).date = datenum(surveyage.textdata(2:end,1),'dd/mm/yyyy');
Y(30).Y    = mean(surveyage.data(:,j(4:5)),2)*100*cj;
Y(30).h    = 0;
Y(30).lag  = 0;
Y(30).age  = 2;
Y(30).hold = 1;

Y(31).type = 'Prevalence 25-, 25-65, 65+ (PHE)'; % Estimated positivity (England)
Y(31).unit = 'percent';
Y(31).U    = 11;
Y(31).date = datenum(surveyage.textdata(2:end,1),'dd/mm/yyyy');
Y(31).Y    = mean(surveyage.data(:,j(6:7)),2)*100*cj;
Y(31).h    = 0;
Y(31).lag  = 0;
Y(31).age  = 3;
Y(31).hold = 0;


% remove NANs, smooth and sort by date
%==========================================================================
M.date  = datestr(min(spm_vec(Y.date)),'dd-mm-yyyy');
M.date  = '01-02-2020';
[Y,YS]  = spm_COVID_Y(Y,M.date,16);

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
pE.sy   = log([1;1;1]);        % coefficients for symptoms
pC.sy   = [1;1;1]/8;           % prior variance(age-specific)

pE.dc   = log([1 1]);          % coefficients for death (28 days)
pC.dc   = [1 1]/8;             % prior variance
pE.ho   = log([1 1]);          % coefficients for admissions
pC.ho   = [1 1]/8;             % prior variance
pE.hc   = log([1 1]);          % coefficients for hospital cases
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
k       = ceil((datenum(date) - datenum(M.date))/64);
pE.tra  = zeros(1,k);          % increases in transmission strength
pC.tra  = ones(1,k)/8;         % prior variance

pE.pcr  = zeros(1,k);          % testing
pC.pcr  = ones(1,k)/8;         % prior variance

pE.mob  = zeros(1,2*k);        % mobility
pC.mob  = ones(1,2*k)/8;       % prior variance

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
M.T       = 28 + datenum(date) - datenum(M.date,'dd-mm-yyyy');
u         = [find(U == 1,1) find(U == 2,1) find(U == 3,1)];
[H,X,~,R] = spm_SARS_gen(Ep,M,[1 2 3]);
spm_SARS_plot(H,X,YS(:,u),[1 2 3])

spm_figure('GetWin','outcomes (1)');
%--------------------------------------------------------------------------
j     = 0;
k     = 0;
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
        if k > 0
            spm_figure('GetWin','outcomes (3)');
        else
            spm_figure('GetWin','outcomes (2)');
            k = k + 1;
        end
        j = 0;
    end
    
end

% time varying parameters
%==========================================================================

% infection fatality ratios (%)
%--------------------------------------------------------------------------
j     = j + 1;
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


%% long-term forecasts Newton(six months from the current data)
%==========================================================================
spm_figure('GetWin','outcomes (4)');

Ep  = DCM.Ep;
Cp  = DCM.Cp;
M   = DCM.M;
M.T = 30*6 + datenum(date) - datenum(M.date,'dd-mm-yyyy');
t   = (1:M.T) + datenum(M.date,'dd-mm-yyyy');

% infection fatality ratios (%)
%--------------------------------------------------------------------------
subplot(2,1,1)
spm_SARS_ci(Ep,Cp,[],19,M); hold on
spm_SARS_ci(Ep,Cp,[],20,M); hold on
ylabel('cases per 100,000'), title('Incidence per 100,000','FontSize',14)
plot(datenum(date)*[1,1],get(gca,'YLim'),':')
legend({'CI per day','actual cases per day','CI per week','confirmed cases per week'})

%% switch windows
%--------------------------------------------------------------------------
spm_figure('GetWin','long-term (1)'); clf

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
[~,~,q] = spm_SARS_ci(Ep,Cp,D,14,M); hold on


% thresholds
%--------------------------------------------------------------------------
% q  = spm_SARS_gen(Ep,M,14); plot(t,q); hold on
%--------------------------------------------------------------------------
u1   = datenum('10-May-2020','dd-mmm-yyyy') - t(1) + 1;
u2   = datenum('10-Aug-2020','dd-mmm-yyyy') - t(1) + 1;
u3   = datenum('31-Aug-2020','dd-mmm-yyyy') - t(1) + 1;
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


%% prevalence and reproduction ratio
%--------------------------------------------------------------------------
spm_figure('GetWin','long-term (2)'); clf

subplot(2,1,1)
i   = find(DCM.U == 4,1);
Rt  = DCM.Y(:,i);
spm_SARS_ci(Ep,Cp,[],11,M); hold on
[~,~,q,c] = spm_SARS_ci(Ep,Cp,Rt,4 ,M); hold on

j   = find(t == datenum(date));
q   = q(j);
d   = sqrt(c{1}(j,j))*1.64;
str = sprintf('Prevalence and reproduction ratio (%s): R = %.2f (CI %.2f to %.2f)',datestr(date,'dd-mmm-yy'),q,q - d,q + d);

% attack rate, herd immunity and herd immunity threshold
%--------------------------------------------------------------------------
[H,~,~,R] = spm_SARS_gen(Ep,M,[4 22 26]);
i         = 16:32;                           % pre-pandemic period
TRN       = [R{1}.Ptrn];                    % transmission risk
R0        = mean(H(i,1));                   % basic reproduction ratio
RT        = R0*TRN(:)/mean(TRN(i));         % effective reproduction ratio
HIT       = 100 * (1 - 1./RT);              % herd immunity threshold
VAC       = H(:,2);                         % percent of people vaccinated
i         = find(H(:,3) > HIT,1);           % date threshold reached
i         = min([i,M.T]);

% Add R0
%--------------------------------------------------------------------------
alpha = datenum('20-Sep-2020','dd-mmm-yyyy');
delta = datenum('20-Mar-2021','dd-mmm-yyyy');
plot(t,RT)
text(alpha,4,'alpha','FontSize',10)
text(delta,4,'delta','FontSize',10)

% add R = 1 and current dateline
%--------------------------------------------------------------------------
plot(get(gca,'XLim'),[1,1],':k')
plot(datenum(date,'dd-mm-yyyy')*[1,1],get(gca,'YLim'),':k')
set(gca,'YLim',[0 8]), ylabel('ratio or percent')
title(str,'FontSize',14)

legend({'CI prevalence','Prevalence (%)','CI R-number','R DCM','R SPI-M','R0'})
legend boxoff
drawnow

% attack rate, herd immunity and herd immunity threshold
%--------------------------------------------------------------------------
subplot(2,1,2)
spm_SARS_ci(Ep,Cp,[],25,M); hold on
spm_SARS_ci(Ep,Cp,[],26,M); hold on

q   = Ep.lnk;
d   = spm_unvec(diag(Cp),Ep);
d   = sqrt(d.lnk)*1.64;
qE  = 100*(1 - exp(q));
qL  = 100*(1 - exp(q + d));
qU  = 100*(1 - exp(q - d));
str = sprintf('Attack rate and immunity: vaccine efficacy %.1f%s (CI %.1f to %.1f)',qE,'%',qL,qU);
q   = Ep.ves;
d   = spm_unvec(diag(Cp),Ep);
d   = sqrt(d.lnk)*1.64;
qE  = 100*(1 - exp(q));
qL  = 100*(1 - exp(q + d));
qU  = 100*(1 - exp(q - d));
sprintf('vaccine efficacy (transmission) %.1f%s (CI %.1f to %.1f)',qE,'%',qL,qU)
leg = sprintf('%s (EIT: %.1f%s)',datestr(t(i),'dd-mmm-yy'),HIT(i),'%');

plot(t,HIT,t,VAC), hold on
plot(t(i)*[1,1],[0,100],':k'), set(gca,'YLim',[0,100])
ylabel('percent'),  title(str,'FontSize',14)
legend({'CI','Attack rate','CI','Population immunity','Effective immunity threshold','Vaccine coverage'})
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

spm_figure('GetWin','outcomes (4)');
savefig(gcf,'Fig4')

spm_figure('GetWin','United Kingdom');
savefig(gcf,'Fig5')

spm_figure('GetWin','long-term (2)');
savefig(gcf,'Fig6')

spm_figure('GetWin','long-term (1)');
savefig(gcf,'Fig7')

% Table
%--------------------------------------------------------------------------
Tab = spm_COVID_table(Ep,Cp,M)

save('DCM_UK.mat','DCM')
cd('C:\Users\karl\Dropbox\Coronavirus')
save('DCM_UK.mat','DCM')

disp('relative transmissibility');
disp(100*TRN(j)/mean(TRN(1:j)))
disp('basic reproduction number');
disp(RT(j))

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

%% Interventions
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(date) - datenum(DCM.M.date,'dd-mm-yyyy');
M.T    = M.T + 180;

u      = 1;
Ep.vef = log(.64);

[Z,X]  = spm_SARS_gen(Ep,M,u,[],3);
spm_SARS_plot(Z,X,S(:,find(U == u(1))),u)


