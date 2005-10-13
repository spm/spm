function t = fieldnames(obj)
% Fieldnames of a file-array object
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: fieldnames.m 253 2005-10-13 15:31:34Z guillaume $

t = {...
    'fname'
    'dim'
    'dtype'
    'offset'
    'scl_slope'
    'scl_inter'
};
