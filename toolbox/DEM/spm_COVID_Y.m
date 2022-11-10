function [Y,S,dates] = spm_COVID_Y(Y,date0,days)
% prepares data array for COVID routines
% FORMAT [Y,S,dates] = spm_COVID_Y(Y,date0)
% Y     - structure array
% date0 - initial date ('dd-mm-yyy')
% days  - number of days over which to average (smooth)
%
% Y     - structure array (time ordered, withough NaNs and smoothed)
% S     - corresponding data matrix
% dates - date numbers from 'dd-mm-yyyy' to last data point
%
%    Y(i).type = datatype (string)
%    Y(i).unit = units (string)
%    Y(i).U    = output index (from spm_SARS_gen);
%    Y(i).date = date number of data points;
%    Y(i).Y    = data points (vector)
%    Y(i).h    = log-precision
%    Y(i).n    = number of data points
%    Y(i).s    = smoothing (days)
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% set up
%==========================================================================
if nargin < 2, date0 = '01-02-2020'; end
if nargin < 3, days  = 7;            end

if ischar(date0)
    date0 = datenum(date0,'dd-mm-yyyy');
end

% check for missing fields
%--------------------------------------------------------------------------
if ~isfield(Y,'lag')
    for i = 1:numel(Y), Y(i).lag = 0; end
end
if ~isfield(Y,'age')
    for i = 1:numel(Y), Y(i).age = 0; end
end
if ~isfield(Y,'h')
    for i = 1:numel(Y), Y(i).h   = 0; end
end

% sort and reorder data
%--------------------------------------------------------------------------
for i = 1:numel(Y)
    
    % remove NaNs
    %----------------------------------------------------------------------
    j         = isfinite(Y(i).Y(:,1));
    Y(i).date = Y(i).date(j);
    Y(i).Y    = Y(i).Y(j,:);
    
    % remove duplicate dates
    %----------------------------------------------------------------------
    [d,j]     = unique(Y(i).date);
    Y(i).date = Y(i).date(j);
    Y(i).Y    = Y(i).Y(j,:);
    
    % remove data prior to initial date
    %----------------------------------------------------------------------
    j         = Y(i).date >= date0;
    Y(i).date = Y(i).date(j);
    Y(i).Y    = Y(i).Y(j,:);
    
    % order data by date
    %----------------------------------------------------------------------
    [d,j]     = sort(Y(i).date);
    Y(i).date = Y(i).date(j);
    Y(i).Y    = Y(i).Y(j,:);
    
end

% smooth data using graph Laplacian (seven day average)
%--------------------------------------------------------------------------
nY    = zeros(1,numel(Y));
for i = 1:numel(Y)
    nY(i)  = numel(Y(i).Y);
    Y(i).n = nY(i);
    Y(i).s = 0;
    
    if Y(i).n
        
        % (AR) precision matrix
        %--------------------------------------------------------------
        Y(i).Q = spm_Q(0,Y(i).n);
        
        % daily timeseries
        %------------------------------------------------------------------
        if mean(diff(Y(i).date)) < 2
            
            % cumulative versus rates
            %--------------------------------------------------------------
            if min(diff(Y(i).Y)) >= 0
                Y(i).Y = cumsum(spm_hist_smooth(gradient(Y(i).Y),days));
            else
                Y(i).Y = spm_hist_smooth(Y(i).Y,days);
            end
            
            % (AR) precision matrix
            %--------------------------------------------------------------
            Y(i).Q = spm_Q(1/2,Y(i).n,1);
            Y(i).s = days;
            
        end
    end
end

% remove low count data
%--------------------------------------------------------------------------
Y     = Y([Y.n] > 1);

% precisions based upon total counts
%--------------------------------------------------------------------------
h     = zeros(1,numel(Y));
for i = 1:numel(Y)
    h(i) = mean(sqrt(Y(i).Y + 1));
end
h     = 4 - log(h);
for i = 1:numel(Y)
    Y(i).h = Y(i).h + h(i);
end

% dates to generate
%--------------------------------------------------------------------------
dates  = date0:max(spm_vec(Y.date));

% data matrix: NaN indicates missing data
%--------------------------------------------------------------------------
S  = NaN(numel(dates),numel(Y));
for i = 1:numel(Y)
    j      = ismember(dates,Y(i).date);
    S(j,i) = Y(i).Y;
end

return




