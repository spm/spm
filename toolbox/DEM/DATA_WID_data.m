function D = DATA_WID_data
% Data retrieval function for COVID modelling
% FORMAT D = DATA_WID_data
%
% n   - number of countries to retain [default: n]
%
% This auxiliary routine retrieves data from comma separated data files
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DATA_WID_data.m 8239 2022-04-09 12:45:02Z karl $


% web options
%--------------------------------------------------------------------------
options = weboptions;
options.Timeout = 20;

% download data and write to CSV files
%--------------------------------------------------------------------------
url = 'https://covid.ourworldindata.org/data/owid-covid-data.csv';
writetable(webread(url,options),'owid-covid-data.csv');

% load data
%--------------------------------------------------------------------------
try
    D  = readtable('owid-covid-data.csv');
catch
    clc, warning('Please load owid-covid-data files into the current working directory')
    return;
end

demo = [ ...
    {'population'                           }
    {'population_density'                   }
    {'median_age'                           }
    {'aged_65_older'                        }
    {'aged_70_older'                        }
    {'gdp_per_capita'                       }
    {'cardiovasc_death_rate'                }
    {'diabetes_prevalence'                  }
    {'hospital_beds_per_thousand'           }
    {'life_expectancy'                      }
    {'human_development_index'              }];

% assemble data structure
%==========================================================================
vname = D.Properties.VariableNames;
State = unique(D{1:end,3});                  % countries
Data  = struct([]);
for k = 1:numel(State)
    j = find(ismember(D{1:end,3},State{k}));
    
    % extract demographic data
    %----------------------------------------------------------------------
    for i = 1:numel(demo)
        Data(k).(demo{i}) = D{j(1),ismember(vname,demo{i})};
    end
    
    
    % extract timeseries
    %----------------------------------------------------------------------
    Data(k).country  = State{k};
    Data(k).date     = D{j,ismember(vname,'date')};
    Data(k).cases    = D{j,ismember(vname,'new_cases')};
    Data(k).death    = D{j,ismember(vname,'new_deaths')};
    Data(k).tests    = D{j,ismember(vname,'new_tests')};
    Data(k).vaccine  = D{j,ismember(vname,'total_vaccinations')};
    Data(k).positive = D{j,ismember(vname,'positive_rate')};
    Data(k).ratio    = D{j,ismember(vname,'reproduction_rate')};
    
    Data(k).vrate    = max(Data(k).vaccine)/Data(k).population;

end



% triage data
%--------------------------------------------------------------------------
i      = zeros(1,numel(Data));
series = {'date','cases','death','vaccine'};
stats  = {'population','aged_70_older'};
for r = 1:numel(Data)
    
    % ensure timeseries are present
    %----------------------------------------------------------------------
    for j = 1:numel(series)
        d = Data(r).(series{j});
        k = datenum(Data(r).date) < datenum('2021-09-01','yyyy-mm-dd');
        if sum(isfinite(d(k))) < 8
            i(r) = 1;
        end
    end
    
    % ensure demographic statistics are present
    %----------------------------------------------------------------------
    for j = 1:numel(stats)
        if isnan(Data(r).(stats{j}))
            i(r) = 1;
        end
    end
end
D     = Data(find(~i));

% remove world
%--------------------------------------------------------------------------
i    = find(ismember({D.country},'World'));
D(i) = [];

if numel(D) ~= 167
    warning('data have changed')
end

return

% Percentage of population currently analysed
%--------------------------------------------------------------------------
load DCM_OWID
pop   = 0;
for r = 1:numel(DCM)
    try, pop = pop + DCM(r).D.population; end
end
i    = find(ismember({Data.country},'World'));
N    = Data(i).population
disp('Perecent popuation')
disp(100 * pop/N)
    

%% order countries using their similarity in terms of new cases and deaths
%==========================================================================
for r = 1:numel(D)
    
    fprintf('%d out of %d\n',r,numel(D));
    disp(D(r).country)
    try, clear Y; end
    
    % create data structure
    %----------------------------------------------------------------------
    Y(1).type = 'New cases';
    Y(1).unit = 'number/day';
    Y(1).U    = 2;
    Y(1).date = datenum(D(r).date);
    Y(1).Y    = D(r).cases;
    Y(1).h    = 0;
    
    Y(2).type = 'Daily deaths';
    Y(2).unit = 'number/day';
    Y(2).U    = 1;
    Y(2).date = datenum(D(r).date);
    Y(2).Y    = D(r).death;
    Y(2).h    = 2;   
    
    % remove NANs, smooth and sort by date
    %----------------------------------------------------------------------
    [Y,S] = spm_COVID_Y(Y,'01-02-2020',16);
    S(isnan(S)) = 0;
    % S       = S/D(r).population;
    
    YY(:,r) = spm_vec(S);
end

[u,s,v] = spm_svd(YY);
[d,i]   = sort(v(:,1));
D       = D(i);

return


