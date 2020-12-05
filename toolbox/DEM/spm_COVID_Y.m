function [Y,S,dates] = spm_COVID_Y(Y,date0,days)
% prepares data array for COVID routines
% FORMAT [Y,S,dates] = spm_COVID_Y(Y,date0)
% Y     - structure array
% date0 - initial date ('dd-mmm-yyy')
% days  - number of days over which to average (smooth)
%
% Y     - structure array (time ordered, withough NaNs and smoothed)
% S     - data matrix
% dates - date numbers from 'dd-mmm-yyyy' to last data point
%
%    Y(i).type = datatype (string)
%    Y(i).unit = units (string)
%    Y(i).U    = output index (from spm_SARS_gen);
%    Y(i).date = date number of data points;
%    Y(i).Y    = data points (vector)
%    Y(i).h    = log-precision
%    Y(i).n    = number of data points
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DEM_COVID_UTLA.m 8005 2020-11-06 19:37:18Z karl $


% set up
%==========================================================================
if nargin < 2, date0 = '01-Feb-2020'; end
if nargin < 3, days  = 7;             end

% remove NANs and sort by date
%--------------------------------------------------------------------------
for i = 1:numel(Y)
    j         = isfinite(Y(i).Y(:,1));
    Y(i).date = Y(i).date(j);
    Y(i).Y    = Y(i).Y(j,:);
    
    [d,j]     = sort(Y(i).date);
    Y(i).date = Y(i).date(j);
    Y(i).Y    = Y(i).Y(j,:);
end

% smooth data using graph Laplacian (seven day average)
%----------------------------------------------------------------------
nY    = zeros(1,numel(Y));
for i = 1:numel(Y)
    nY(i)  = numel(Y(i).Y);
    Y(i).n = nY(i);
    if max(diff(Y(i).date)) < 2
        Y(i).Y = spm_hist_smooth(Y(i).Y,days);
    end
end

% precisions based upon total counts
%----------------------------------------------------------------------
h     = zeros(1,numel(Y));
for i = 1:numel(Y)
    h(i) = sum(Y(i).Y);
end
h   = log(h);
h   = mean(h) - h;
for i = 1:numel(Y)
    Y(i).h = Y(i).h + h(i);
end

% dates to generate
%--------------------------------------------------------------------------
dates  = datenum(date0,'dd-mmm-yyyy'):max(spm_vec(Y.date));

% data matrix (smooth): NaN indicates missing data
%----------------------------------------------------------------------
S  = NaN(numel(dates),numel(Y));
for i = 1:numel(Y)
    j      = ismember(dates,Y(i).date);
    S(j,i) = Y(i).Y;
end




