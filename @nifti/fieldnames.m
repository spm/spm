function t = fieldnames(obj)
% Fieldnames of a NIFTI-1 object
% _______________________________________________________________________
% %W% John Ashburner %E%

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
