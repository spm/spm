function res = nfrequencies(this)
% Method for getting the number of frequencies for TF data
% FORMAT res = nsamples(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: nfrequencies.m 1270 2008-03-28 14:35:16Z stefan $

if ~strcmp(transformtype(this), 'TF')
    res = [];
else
    res = length(this.transform.frequencies);
end