function Data = DATA_WID_data
% Data retrieval function for COVID modelling
% FORMAT data = DATA_WID_data
%
% n   - number of countries to retain [default: n]
%
% This auxiliary routine retrieves data from comma separated data files
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DATA_COVID_JHU.m 8129 2021-08-02 18:08:36Z karl $


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
    
    % check
    %----------------------------------------------------------------------
    Data(k).country  = State{k};
    Data(k).date     = D{j,ismember(vname,'date')};
    Data(k).cases    = D{j,ismember(vname,'new_cases')};
    Data(k).death    = D{j,ismember(vname,'new_deaths')};
    Data(k).tests    = D{j,ismember(vname,'new_tests')};
    Data(k).vaccine  = D{j,ismember(vname,'people_vaccinated')};
    Data(k).positive = D{j,ismember(vname,'positive_rate')};
    Data(k).ratio    = D{j,ismember(vname,'reproduction_rate')};

    % check
    %----------------------------------------------------------------------
    for i = 1:numel(demo)
        Data(k).(demo{i}) = D{j(1),ismember(vname,demo{i})};
    end

end

return


