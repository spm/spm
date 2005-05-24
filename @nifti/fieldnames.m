function t = fieldnames(obj)
% Fieldnames of a NIFTI-1 object
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: fieldnames.m 174 2005-05-24 11:03:32Z john $


if isfield(obj.hdr,'magic')
    t = {...
        'dat'
        'mat'
        'mat_intent'
        'mat0'
        'mat0_intent'
        'intent'
        'diminfo'
        'timing'
        'descrip'
        'cal'
        'aux_file'
    };
else
    error('This should not happen.');
end;
