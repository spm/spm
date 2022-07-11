function res = size(this, dim)
% returns the dimensions of the data matrix
% FORMAT res = size(this, dim))
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if ~strncmpi(transformtype(this), 'TF', 2)
    res = [nchannels(this), nsamples(this), ntrials(this)];
else
    res = [nchannels(this), nfrequencies(this), nsamples(this), ntrials(this)];
end

if nargin > 1
    res = res(dim);
end
