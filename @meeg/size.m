function res = size(this, dim)
% returns the dimensions of the data matrix
% FORMAT res = size(this, dim))
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: size.m 5078 2012-11-25 15:08:05Z vladimir $


if ~strncmpi(transformtype(this), 'TF', 2)
    res = [nchannels(this), nsamples(this), ntrials(this)];
else
    res = [nchannels(this), nfrequencies(this), nsamples(this), ntrials(this)];
end

if nargin > 1
    res = res(dim);
end