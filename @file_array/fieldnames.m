function t = fieldnames(obj)
% Fieldnames of a file-array object
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: fieldnames.m 7147 2017-08-03 14:07:01Z spm $


t = {...
    'fname'
    'dim'
    'dtype'
    'offset'
    'scl_slope'
    'scl_inter'
};
