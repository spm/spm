function spm_eeg_inv_image_display(varargin)
% Display an interpolated 3D image or mesh of a contrast or window
%
% FORMAT D = spm_eeg_inv_image_display(D,val)
% Input:
% D        - input data struct (optional)
%__________________________________________________________________________
% Copyright (C) 2007-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_image_display.m 6405 2015-04-14 15:13:02Z guillaume $


%-Checks
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

%-EEG and underlay files
%--------------------------------------------------------------------------
try
    sMRI = D.inv{val}.mesh.wmMRI; spm_vol(sMRI);
catch
    sMRI = fullfile(spm('dir'),'canonical','single_subj_T1.nii');
end
wEEG = D.inv{val}.contrast.fname{D.con};

%-Display
%--------------------------------------------------------------------------
spm_figure('Clear','Graphics');
switch lower(char(D.inv{val}.contrast.format))
    case 'image'
        spm_check_registration(sMRI);
        spm_orthviews('addcolouredimage',1,[wEEG, ',1'],[1 0 0]);
        spm_orthviews('addcolourbar',1,1);
        spm_orthviews('Redraw');
    case 'mesh'
        ax = subplot(2,1,1,'parent',spm_figure('GetWin','Graphics'));
        spm_mesh_render('Disp',wEEG,'parent',ax);
    otherwise
        error('Unknown data format.');
end
