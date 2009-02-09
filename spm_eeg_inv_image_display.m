function spm_eeg_inv_image_display(varargin)
% spm_orthviews of an interpolated 3D image of a contrast or window
%
% FORMAT D = spm_eeg_inv_image_display(D,val)
% Input:
% D        - input data struct (optional)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout, Stefan Kiebel & Karl Friston
% $Id: spm_eeg_inv_image_display.m 2720 2009-02-09 19:50:46Z vladimir $


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
spm_check_registration(sMRI);
spm_orthviews('addcolouredimage',1,wEEG,[1 0 0]);
spm_orthviews('addcolourbar',1,1);
spm_orthviews('Redraw');
rotate3d off;
