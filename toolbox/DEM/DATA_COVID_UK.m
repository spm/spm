function [Y,R] = DATA_COVID_UK(country)
% Data retrieval function for COVID modelling (UK)
% FORMAT [Y,R] = DATA_COVID_UK(country)
%
% Y  - daily deaths and confirmed cases
% R  - daily test rates for the UK
%
% country - optional country
%
% This auxiliary routine retrieves data from comma separated data files
% that can be downloaded from:
% https://github.com/CSSEGISandData/COVID-19/
% https://github.com/tomwhite/covid-19-uk-data
%
%     time_series_covid19_confirmed_global.csv
%     time_series_covid19_deaths_global.csv
%     time_series_covid19_recovered_global.csv
%     covid-19-tests-uk
%
% It augments these data with population sizes from the United Nations,
% returning the following data structure:
%
% Data(k).country - country
% Data(k).pop     - population size
% Data(k).lat     - latitude
% Data(k).long    - longitude
% Data(k).date    - date when more than one case was reported
% Data(k).cases   - number of cases,  from eight days prior to first cases
% Data(k).death   - number of deaths, from eight days prior to first cases
% Data(k).recov   - number recovered, from eight days prior to first cases
% Data(k).days    - number of days in timeseries
% Data(k).cum     - cumulative number of deaths
%
% Population data from (cite as):
% United Nations, Department of Economic and Social Affairs, Population
% Division (2019). World Population Prospects 2019, Online Edition. Rev. 1.
%
% Please see the main body of the script for a description of the graphical
% outputs provided when the routine is called with at an output argument.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% https://coronavirus.data.gov.uk/coronavirus-cases_latest.csv
% https://coronavirus.data.gov.uk/coronavirus-deaths_latest.csv


% download data if required
%==========================================================================
url = 'https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/';
urlwrite([url,'time_series_covid19_confirmed_global.csv'],'time_series_covid19_confirmed_global.csv');
urlwrite([url,'time_series_covid19_deaths_global.csv'],   'time_series_covid19_deaths_global.csv');
urlwrite([url,'time_series_covid19_recovered_global.csv'],'time_series_covid19_recovered_global.csv');

url = 'https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/';
urlwrite([url,'covid-19-totals-uk.csv'],'covid-19-totals-uk.csv');

% defaults
%--------------------------------------------------------------------------
if nargin < 1, country = 'United Kingdom'; end


% load data from https://github.com/CSSEGISandData/COVID-19/
%--------------------------------------------------------------------------
try
    C  = importdata('time_series_covid19_confirmed_global.csv');
    D  = importdata('time_series_covid19_deaths_global.csv'   );
    R  = importdata('time_series_covid19_recovered_global.csv');
    T  = importdata('covid-19-totals-uk.csv');
catch
    clc, warning('Please load csv files into the current working directory')
    help DATA_COVID_UK
    return;
end

% dates of JHU data
%--------------------------------------------------------------------------
date = D.textdata(1,5:end);                     % date

% assemble data structure
%==========================================================================
s   = 7;                                          % data smoothing (days)

% confirmed cases
%--------------------------------------------------------------------------
Ci  = logical(ismember(C.textdata(2:end,2),country));
CY  = sum(C.data(Ci,3:end),1)';

% confirmed deaths
%--------------------------------------------------------------------------
Di  = logical(ismember(D.textdata(2:end,2),country));
DY  = sum(D.data(Di,3:end),1)';

% recovered
%--------------------------------------------------------------------------
Ri  = logical(ismember(R.textdata(2:end,2),country));
RY  = sum(R.data(Ri,3:end),1)';

% Total tests, and align first date
%--------------------------------------------------------------------------
TY  = T.data(:,1);
TT  = T.textdata(2:end,1);             disp(TT{1})
Ti  = find(ismember(date,'1/25/20'));  disp(date{Ti})

% take gradients of cumulative rates
%--------------------------------------------------------------------------
DY  = gradient(spm_conv(DY(Ti:end),s));
CY  = gradient(spm_conv(CY(Ti:end),s));
RY  = gradient(spm_conv(RY(Ti:end),s));
TY  = gradient(spm_conv(TY(1 :end),s));
TY(isnan(TY)) = max(TY);

% assemble response variables
%==========================================================================
Y   = [DY CY];
R   = TY;





