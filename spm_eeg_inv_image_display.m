function spm_eeg_inv_image_display(varargin)

%======================================================================
% spm_orthviews of an interpolated 3D image of a contrast or winow
%
% FORMAT D = spm_eeg_inv_image_display(D,val)
% Input:
% D		   - input data struct (optional)
%
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout, Stefan Kiebel & Karl Friston
% $Id: spm_eeg_inv_image_display.m 459 2006-02-24 11:55:50Z jeremie $


% checks
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% display
%==========================================================================
Fgraph  = spm_figure('GetWin','Graphics'); clf(Fgraph), figure(Fgraph)

% EEG and sMRI files
%--------------------------------------------------------------------------
try
    sMRI = D.inv{val}.mesh.wmMRI;
    wEEG = D.inv{val}.contrast.fname;
catch
    sMRI = fullfile(spm('dir'),'canonical','single_subj_T1.nii');
    wEEG = D.inv{val}.contrast.fname;
end

job{1}.util{1}.checkreg.data{1} = wEEG;
job{1}.util{1}.checkreg.data{2} = sMRI;
spm_jobman('run',job);

