function t = fieldnames(obj)
% Fieldnames of a file-array object
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: fieldnames.m 174 2005-05-24 11:03:32Z john $

t = {...
    'fname'
    'dim'
    'dtype'
    'offset'
    'scl_slope'
    'scl_inter'
};
