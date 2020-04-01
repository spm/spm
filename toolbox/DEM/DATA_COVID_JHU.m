function data = DATA_COVID_JHU
% Data retrieval function for COVID modelling
% FORMAT data = DATA_COVID_JHU
%
% This auxiliary routine retrieves data from comma separated data files
% that can be downloaded from:
% https://github.com/CSSEGISandData/COVID-19/
%
%     time_series_covid19_confirmed_global.csv
%     time_series_covid19_deaths_global.csv
%     time_series_covid19_recovered_global.csv
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
% $Id: DATA_COVID_JHU.m 7810 2020-04-01 13:58:56Z spm $


% load data from https://github.com/CSSEGISandData/COVID-19/
%--------------------------------------------------------------------------
try
    C  = importdata('time_series_covid19_confirmed_global.csv');
    D  = importdata('time_series_covid19_deaths_global.csv'   );
    R  = importdata('time_series_covid19_recovered_global.csv');
catch
    clc, warning('Please load csv files into the current working directory')
    help DATA_COVID_JHU
end

N  = importdata('population.xlsx');          % population size


% preliminary extraction
%--------------------------------------------------------------------------
date       = D.textdata(1,5:end);            % date
Location   = D.textdata(2:end,2);            % location by country
State      = unique(Location);               % common countries
Npop       = N.data;
Country    = N.textdata;

% ensure consistency between covid timeseries and population data
%--------------------------------------------------------------------------
i          = logical(ismember(Country,{'United States of America','US'}));
Country{i} = 'US';


% assemble data structure
%==========================================================================
Data  = struct([]);
k     = 1;
for i = 1:numel(State)
    j = find(ismember(Country,State{i}));
    if numel(j) 
        
        % confirmed cases
        %------------------------------------------------------------------
        Ci  = logical(ismember(C.textdata(2:end,2),State{i}));
        CY  = sum(C.data(Ci,3:end),1)';

        % confirmed deaths
        %------------------------------------------------------------------
        Di  = logical(ismember(D.textdata(2:end,2),State{i}));
        DY  = sum(D.data(Di,3:end),1)';
        
        % recovered
        %------------------------------------------------------------------
        Ri  = logical(ismember(R.textdata(2:end,2),State{i}));
        RY  = sum(R.data(Ri,3:end),1)';
        
        % manual check
        %------------------------------------------------------------------
        % if ismember(State{i},'China'), keyboard, end
        
        % save from first reported case
        %------------------------------------------------------------------
        d   = find(cumsum(CY) > 1,1);
        l   = find(ismember(C.textdata(2:end,2),State{i}),1);
        
        if sum(CY) > 32
            Data(k).country = State{i};
            Data(k).pop     = Npop(j)*1e3;
            Data(k).lat     = C.data(l,1);
            Data(k).long    = C.data(l,2);
            Data(k).date    = date{d};
            Data(k).cases   = gradient(spm_conv([zeros(8,1); CY(d:end)],2));
            Data(k).death   = gradient(spm_conv([zeros(8,1); DY(d:end)],2));
            Data(k).recov   = gradient(spm_conv([zeros(8,1); RY(d:end)],2));
            Data(k).days    = numel(Data(k).cases);
            Data(k).cum     = sum(Data(k).death);
            k = k + 1;
        end
        
    end
end

% Illustrate some features of the data, focusing on early cases
%==========================================================================
% This figure provides a brief survey of the timeseries used for subsequent
% modelling, with a focus on the early trajectories of mortality. The upper
% left panel shows the distribution, over countries, of the number of days
% after more than one case was reported. At the time of writing, a
% substantial number of countries witnessed an outbreak lasting for more
% than 60 days. The upper left panel plots the total number of deaths
% against the durations in the left panel. Those countries whose outbreak
% started earlier have greater cumulative deaths; however, within this
% group, there is no clear correlation between population size and total.
% The middle left panel plots the new deaths reported (per day) over a
% 48-day period following the first report of more than one case. The
% colours of the lines denote different countries. These countries are
% listed in the lower left panel, which plots the cumulative death rate.
% China is clearly the first country to be severely affected, with
% remaining countries evincing an accumulation of deaths some 30 days after
% China. Interestingly, there is little correlation between the total
% number of deaths and population size. The middle right panel is a
% logarithmic plot of the total deaths against population size in the
% initial (48-day) period. However, there is a stronger correlation between
% the total number of cases reported (within the first 48 days) and the
% total number of deaths.

% duration of pandemic
%--------------------------------------------------------------------------
T   = spm_cat({Data.days});
N   = spm_cat({Data.cum});

% retain countries with over 48 days and 16 deaths
%--------------------------------------------------------------------------
t    = 48;
i    = logical(T > t & N > 16);
data = Data(i);
for i = 1:numel(data)
    death(:,i) = data(i).death(1:t);
    cases(:,i) = data(i).cases(1:t);
end
[d,j] = sort(sum(death),'descend');
death = death(:,j);
cases = cases(:,j);
data  = data(j);
t     = 1:t;

% plot sorting, unless there is an outcome argument specified
%--------------------------------------------------------------------------
if nargout, return, end
    
spm_figure('GetWin','Data');

subplot(3,2,1), hist(T,32,'Edgecolor','none')
title('Duration (days)','Fontsize',16)
xlabel('days'),ylabel('number of countries'), axis square

subplot(3,2,2), plot(T,N,'.','MarkerSize',32,'Color',[0.5 0.5 1])
title('Total deaths','Fontsize',16)
xlabel('days'),ylabel('total deaths'), axis square

% compare initial statistics
%--------------------------------------------------------------------------
subplot(3,2,3), plot(t,death)
title('Death rate','Fontsize',16)
xlabel('days'),ylabel('new cases')
axis square

pop = [data.pop];
subplot(3,2,4), loglog(pop,sum(death),'.','MarkerSize',32,'Color',[0.8 0.8 1])
title('Population and deaths','Fontsize',16)
xlabel('population'),ylabel('total deaths (within 48 days)')
axis square

subplot(3,2,5), plot(t,cumsum(death,1))
title('Cumulative deaths','Fontsize',16)
xlabel('days'),ylabel('new cases')
axis square, legend({data(1:14).country}), legend('boxoff')

subplot(3,2,6), loglog(sum(death),sum(cases),'.','MarkerSize',32,'Color',[0.8 0.8 1])
title('Cases and deaths','Fontsize',16)
xlabel('cumulative deaths'),ylabel('total cases (within 48 days)')
axis square

% label countries
%--------------------------------------------------------------------------
for i = 1:numel(data)
    subplot(3,2,6), text(sum(death(:,i)),sum(cases(:,i)),data(i).country,'FontSize',9)
    subplot(3,2,4), text(pop(i),sum(death(:,i)),         data(i).country,'FontSize',9)
end
