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
% $Id: spm_eeg_inv_image_display.m 1143 2008-02-07 19:33:33Z spm $


% checks
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% display
%==========================================================================
Fgraph  = spm_figure('GetWin','Graphics'); clf(Fgraph), figure(Fgraph)

% EEG and sMRI files
%--------------------------------------------------------------------------
try
    sMRI = D.inv{val}.mesh.wmMRI; spm_vol(sMRI);
    wEEG = D.inv{val}.contrast.fname{D.con};
catch
    sMRI = fullfile(spm('dir'),'canonical','single_subj_T1.nii');
    wEEG = D.inv{val}.contrast.fname{D.con};
end

job{1}.util{1}.checkreg.data{1} = wEEG;
job{1}.util{1}.checkreg.data{2} = sMRI;
spm_jobman('run',job);


