function res = size(this, varargin)
% returns the dimensions of the data matrix
% FORMAT res = size(this, dim))
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: size.m 5061 2012-11-16 11:15:50Z vladimir $


if ~strncmpi(transformtype(this), 'TF', 2)
    res = [nchannels(this), nsamples(this), ntrials(this)];
else
    res = [nchannels(this), nfrequencies(this), nsamples(this), ntrials(this)];
end

