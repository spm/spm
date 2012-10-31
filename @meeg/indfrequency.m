function res = indfrequency(this, f)
% Method for getting the index closest to given frequency
% FORMAT  res = indfrequency(this, f)
% this       - MEEG object
% f          - vector of frequencies (in Hz)
%
% res        - vector of sample indices matching indices
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: indfrequency.m 5025 2012-10-31 14:44:13Z vladimir $

if ~strncmpi(transformtype(this), 'TF',2)
    error('Only TF datasets are supported');
end

res = NaN(1,length(f));
fdiff = mean(diff(frequencies(this)));
if nsamples(this) > 0
    F = frequencies(this);
    for i = 1:length(f)
        [m, res(i)] = min(abs(F-f(i)));
        if m > fdiff
            warning('Could not find an index matching the requested frequency %d Hz', f(i));
            res(i) = NaN;
        end
    end
end
