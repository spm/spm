function [vbmix logev mix] = spm_mix_demo1d (data, maxcomps, verbosity)
% Demonstrate use of spm_mix on 1D data
% FORMAT [vbmix logev mixdata] = spm_mix_demo1d (data, maxcomps, plotfits)
%
% data      - either scalar number of clusters to simulate or your own data
% maxcomps  - maximum number of components in mixture model to consider
% verbosity - 0 = silent, 1 = basic output (with figures), 2 = full output
%
% vbmix     - cell array of fitted mixtures for all numbers of components
% logev     - log evidence for each number of components
% mix       - mix structure for simulated mixtures if scalar data given
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Will Penny & Ged Ridgway
% $Id: spm_mix_demo1d.m 3997 2010-07-15 12:38:24Z ged $

if nargin < 3
    verbosity = 1;
end
if nargin < 2
    maxcomps = 5;
end
if nargin < 1
    data = 3;
end

if isscalar(data)
    mix.m = data;
    if verbosity > 0
        fprintf('Simulating data with %d clusters\n', mix.m)
    end
    means = 0:5:5*(mix.m - 1);
    for m = 1:mix.m
        mix.state(m).prior = 1 / mix.m;
        mix.state(m).m = means(m);
        mix.state(m).C = 1;
    end
    N = 50 * mix.m; % Number of data points
    data = spm_samp_mix(mix, N);
else
    data = data(isfinite(data));
    mix = 'Own data given';
end

logev = nan(maxcomps, 1);
vbmix = cell(maxcomps, 1);
for m = 1:maxcomps
    if verbosity > 0
        fprintf('Fitting mixture model with %d components\n', m);
    end
    vbmix{m} = spm_mix(data, m, verbosity > 1);
    if verbosity > 0
        figure
        spm_mix_plot1d (data, vbmix{m})
        title('Fitted model')
    end
    logev(m) = vbmix{m}.fm;
end
logev = logev-min(logev);

if verbosity > 0
    if isstruct(mix)
        figure;
        spm_mix_plot1d (data, mix);
        title('True generative model')
    end
    figure
    bar(logev);
    xlabel('Number of mixture components');
    ylabel('Log Evidence');
end

if nargout == 0 && verbosity > 0 % assume wanted plots without other output
    clear vbmix
end
