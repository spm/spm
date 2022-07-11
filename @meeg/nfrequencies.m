function res = nfrequencies(this)
% Method for getting the number of frequencies for TF data
% FORMAT res = nsamples(this)
%__________________________________________________________________________

% Stefan Kiebel
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if ~strncmp(transformtype(this), 'TF',2)
    res = [];
else
    res = length(this.transform.frequencies);
end
