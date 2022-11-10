function [Y,Data] = DATA_COVID_US
% Data retrieval function for COVID modelling
% FORMAT [Y,Data] = DATA_COVID_US
%
% Y(:,i,1) = Data(i).death;
% Y(:,i,2) = Data(i).cases;
%
% This auxiliary routine retrieves data from comma separated data files
% that can be downloaded from:
% https://github.com/CSSEGISandData/COVID-19/
%
%     time_series_covid19_confirmed_US.csv
%     time_series_covid19_deaths_US.csv
%
% Data(k).state   - State
% Data(k).pop     - population size
% Data(k).lat     - latitude
% Data(k).long    - longitude
% Data(k).date    - date of 8th day
% Data(k).cases   - number of cases
% Data(k).death   - number of deaths
% Data(k).cum     - cumulative number of deaths
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% download data
%==========================================================================
url = 'https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/';
urlwrite([url,'time_series_covid19_confirmed_US.csv'],'time_series_covid19_confirmed_US.csv');
urlwrite([url,'time_series_covid19_deaths_US.csv'],   'time_series_covid19_deaths_US.csv');

% load data from https://github.com/CSSEGISandData/COVID-19/
%--------------------------------------------------------------------------
try
    C  = importdata('time_series_covid19_confirmed_US.csv',',',1);
    D  = importdata('time_series_covid19_deaths_US.csv',',',1   );
catch
    clc, warning('Please load csv files into the current working directory')
    help DATA_COVID_US
end


% preliminary extraction
%--------------------------------------------------------------------------
Region = D.textdata(2:end,7);            % location by country
CD     = D.data(:,2:end);
CC     = C.data(:,1:end);
Pop    = D.data(:,1);

for i = 1:length(Pop)
    try
        Lat(i) = str2num(D.textdata{i + 1,9});
        Lon(i) = str2num(D.textdata{i + 1,10});
    catch
        Lat(i) = 0;
        Lon(i) = 0;
    end
end
State  = unique(Region);                 % common States


% assemble data structure
%==========================================================================
Data  = struct([]);
s     = 2;                                          % data smoothing (days)
T     = 16;                                         % days to skip
k     = 1;
for i = 1:numel(State)
    j = find(ismember(Region,State{i}));
    
    % confirmed cases
    %----------------------------------------------------------------------
    CY  = sum(CC(j,:),1)';
    
    % confirmed deaths
    %----------------------------------------------------------------------
    DY  = sum(CD(j,:),1)';
    
    % population
    %----------------------------------------------------------------------
    PY  = sum(Pop(j));
    
    % latitude and longitude
    %----------------------------------------------------------------------
    lat = mean(Lat(j));
    lon = mean(Lon(j));
    
    % create data structure
    %----------------------------------------------------------------------
    if DY(end) > 32
        Data(k).state  = State{i};
        Data(k).pop    = PY;
        Data(k).lat    = lat;
        Data(k).lon    = lon;
        Data(k).date   = '20/1/2020';
        Data(k).cases  = gradient(spm_conv(CY(T:end),s));
        Data(k).death  = gradient(spm_conv(DY(T:end),s));
        Data(k).cum    = sum(Data(k).death);
        Data(k).first  = find(CY,1,'first');
        k = k + 1;
    end
    
end

% rank and compute geographical distances
%==========================================================================
[d,i] = sort([Data(:).cum],'descend');
Data  = Data(i);
for i = 1:numel(Data)
    
    % geographical distances
    %----------------------------------------------------------------------
    Data(i).dis = ([Data.lat]' - Data(i).lat).^2 + ...
                  ([Data.lon]' - Data(i).lon).^2;
              
    % onset distances
    %----------------------------------------------------------------------
    Data(i).ons = ([Data.first]' - Data(i).first).^2;
    
    % data
    %----------------------------------------------------------------------
    Y(:,1,i)    = Data(i).death;
    Y(:,2,i)    = Data(i).cases;
    
end

% remove negative values
%--------------------------------------------------------------------------
Y(Y < 0) = 0;
