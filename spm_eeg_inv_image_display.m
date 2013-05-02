function spm_eeg_inv_image_display(varargin)
% Display an interpolated 3D image of a contrast or window
%
% FORMAT D = spm_eeg_inv_image_display(D,val)
% Input:
% D        - input data struct (optional)
%__________________________________________________________________________
% Copyright (C) 2007-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_image_display.m 5461 2013-05-02 19:01:57Z vladimir $


% checks
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% EEG and sMRI files
%--------------------------------------------------------------------------
try
    sMRI = D.inv{val}.mesh.wmMRI; spm_vol(sMRI);
catch
    sMRI = fullfile(spm('dir'),'canonical','single_subj_T1.nii');
end
wEEG = D.inv{val}.contrast.fname{D.con};

% display
%--------------------------------------------------------------------------
spm_figure('Clear','Graphics');
spm_check_registration(sMRI);
spm_orthviews('addcolouredimage',1,[wEEG, ',1'],[1 0 0]);
spm_orthviews('addcolourbar',1,1);
spm_orthviews('Redraw');
