function [X] = DEM_COVID_MASKS(i)
% Returns the probability of mask wearing at a particular day
% FORMAT [X] = DEM_COVID_MASKS(i)
% = datenum('dd-mmm-yyyy')
%__________________________________________________________________________


% Percentage of adults that have used a face covering when outside their home in the past seven days	
%---------------------------------------------------------------------------------------------------
x = [28 28 30 39 44 43 52 61 71 84 96 96 95 96 95 96 97 98 98 96 98 97 98 97 97 97 97 97 97 97 96 96 95 94 95 96 96 95 96 97 97 97 96 97 97 97 98 97 97 97 97 96 96 97 97 96 94 95 95 92 90 89 90 89 88 86 82 83 85 84 94 96 95 95 88 83 74 68 63 65 57 54 48 38 35 36 34 35 30 26 25 29 29 24 21 21 22 22 21 21 19 18 16 ];
date0 = datenum('07-Jun-2020');
date1 = datenum('22-Mar-2023');

% interpolate ONS data
%--------------------------------------------------------------------------
t = fix(linspace(date0,date1,numel(x)));
T = date0:date1;
X = resample(x/100,numel(T),numel(t));

% return prevalence of mask wearing
%--------------------------------------------------------------------------
if i < T(1)
    X = 0;
elseif i > (T(end) - 4)
    X = X(end - 4);
else
    X = X(ismember(T,i));
end
