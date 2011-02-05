function [] = spm_mix_plot1d (data, mix, rng, nPoints)
% Plot component densities and mixture density for 1D mixture model
% FORMAT [] = spm_mix_plot1d (data, mix, rng, nPoints)
%
% data     - Optional original data, from which histogram is plotted, or []
% mix      - Mixture model data structure from spm_mix
% rng      - [xmin xmax], defaulting to [min(data) max(data)] if data given
% nPoints  - Number of points covering rng; default=100
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Will Penny & Ged Ridgway
% $Id: spm_mix_plot1d.m 4195 2011-02-05 18:39:13Z ged $

if nargin < 4 || isempty(nPoints), nPoints = 100; end
if nargin < 3 || isempty(rng)
    if ~isempty(data)
        rng = [min(data) max(data)];
    else
        error('Please either specify range or give data to derive it from')
    end
end

x = linspace(min(rng), max(rng), nPoints);
M = mix.m;
pat = sprintf('Comp. %%0%dd', 1+floor(log10(M)));
pdfs = zeros(M + 1, nPoints);
legs = cell(M+1, 1);
if M == 1
    mix.state(1).prior = 1;
end
for m = 1:M
    pdfs(m, :) = mix.state(m).prior * ...
        spm_Npdf(x, mix.state(m).m, mix.state(m).C);
    legs{m} = sprintf(pat, m);
end
pdfs(M+1, :) = sum(pdfs);
legs{M+1} = 'Mixture PDF';

if ~isempty(data)
    subplot(2,1,1)
end
plh = plot(x, pdfs);
mixcol = get(plh(end), 'Color');
set(plh(end), 'LineStyle', '--')
legend(legs, 'Location', 'BestOutside');

if ~isempty(data)
    xl = xlim;
    subplot(2,1,2)
    [c b] = hist(data, round(sqrt(numel(data))));
    bar(b, c / sum(c) / mean(diff(b)), 'k');
    if exist('ksdensity', 'file')
        hold on;
        f = ksdensity(data, x);
        plot(x, f, '--', 'Color', mixcol);
        legend('Histogram', 'Parzen PDF', 'Location', 'BestOutside');
    end
    xlim(xl);
    subplot(2,1,1); % (so any subsequent title command goes over this plot)
end
