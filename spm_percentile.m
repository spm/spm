function [y] = spm_percentile(data, p)
% Compute one or more percentiles from data
% FUNCTION [y] = spm_percentile(data, p)
% data - arbirarily sized input data (from which NaNs will be excluded)
% p    - scalar or n-vector of percentage values (from 0 to 100)
%        if not specified, p defaults to all quartiles: [0 25 50 75 100]
%
% y    - scalar or n-vector of corresponding percentiles
%
% Note that percentiles are computed over all data, not along the first or
% specified dimension (unlike prctile from the MATLAB Statistics Toolbox).
%
% Example:
%  spm_summarise(vols, 'all', @spm_percentile) % quartiles of images
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: spm_percentile.m 4013 2010-07-22 17:12:45Z guillaume $

% The algorithm used is that described by NIST, with x = k+d = 1+p(N-1)/100
% http://www.itl.nist.gov/div898/handbook/prc/section2/prc252.htm
%
% This choice apparently matches Excel, but not MATLAB's prctile, though 
% the differences are typically small (and are zero for min, median, max).
%
% This algorithm was chosen because it requires no special handling of 0 or
% 100, and lets 0 and 1 percentiles differ even with less than 100 samples.
% It also has the appealing property of returning uniformly spaced values
% for uniformly spaced percentiles of uniformly spaced data. For example:
%  x = 1:11;
%  p = 0:25:100;
%  diff([prctile(x, p(:)) spm_percentile(x, p) spm_percentile_nist(x, p)])
% gives constant differences only for spm_percentile.


if nargin < 2
    p = 0:25:100;
end
if any(p < 0 | p > 100)
    error('Percentage values outside 0 <= p <= 100 are not allowed')
end
if all(p < 1) && any(p > 0)
    warning('SPM:invalidPercentageValues', ...
        'Percentage values are below 1, but values up to 100 are expected')
end

p    = p(:);
data = data(:);

data        = sort(data(~isnan(data)));
N           = length(data);
data(end+1) = data(end); % to handle special case of 100th percentile below

x = 1 + p * (N - 1) / 100;
k = floor(x);
d = x - k;

y = (1 - d) .* data(k) + d .* data(k + 1);
