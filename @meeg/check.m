function [ok, this] = check(this, option)
% Method that performs integrity checks of the meeg object
% and its readiness for particular purposes.
% FORMAT  this = check(this, option)
% IN
% option - 'basic' (default) - just check the essential fields
%          'sensfid' - also checks sensor and fiducial definitions
% OUT
% ok - 1 - OK, 0- failed
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: check.m 2720 2009-02-09 19:50:46Z vladimir $

[ok, this] = checkmeeg(struct(this), option);


this = meeg(this);

if ok && strcmp(option, '3d')
    if ~ismember(modality(this), {'EEG', 'MEG', 'Multimodal'})
        warning('Unsupported modality for 3D source reconstruction');
        ok = 0;
    end
end

if ok && strcmp(option, 'dcm')
    if ~ismember(modality(this, 0), {'EEG', 'MEG', 'MEGPLANAR', 'Multimodal', 'LFP'})
        warning('Unsupported modality for DCM');
        ok = 0;
    end
end