function t = fieldnames(obj)
% Fieldnames of a NIFTI-1 object
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: fieldnames.m 253 2005-10-13 15:31:34Z guillaume $


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
