% This script defines spatial confounds and adds them to MEEG dataset. This
% This functionality will be further developed in the future. For now spatial
% confounds are loaded from a BESA *.bsa file.
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_spatial_confounds.m 2333 2008-10-13 13:19:22Z vladimir $

D = spm_eeg_load;

fbsa = spm_select(1, '\.bsa$', 'Select BESA *.bsa file');

D = sconfounds(D, spm_eeg_read_bsa(fbsa));

save(D);