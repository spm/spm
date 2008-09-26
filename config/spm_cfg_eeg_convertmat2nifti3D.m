function S = spm_cfg_eeg_convertmat2nifti3D
% configuration file for writing voxel-based images from SPM M/EEG format,
% as a time-series of 2Dimages
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_convertmat2nifti3D.m 2201 2008-09-26 10:14:50Z stefan $

Fname = cfg_files;
Fname.tag = 'Fname';
Fname.name = 'File Names';
Fname.filter = 'mat';
Fname.num = [1 inf];
Fname.help = {'Select the M/EEG mat file.'};

n = cfg_entry;
n.tag = 'n';
n.name = 'Output dimension';
n.strtype = 'r';
n.num = [1 1];
n.val = {64};
n.help = {'Enter the Output image dimension'};

pixsize = cfg_entry;
pixsize.tag = 'pixsize';
pixsize.name = 'Voxel size';
pixsize.strtype = 'r';
pixsize.num = [1 1];
pixsize.val = {3};
pixsize.help = {'Enter the approximate voxel output size'};

yes = cfg_const;
yes.tag = 'yes';
yes.name = 'Interpolate bad channels';
yes.val = {1};

no = cfg_const;
no.tag = 'no';
no.name = 'Mask out bad channels';
no.val = {1};

Interpolate = cfg_choice;
Interpolate.tag = 'interpolate_bad';
Interpolate.name = 'Interpolate';
Interpolate.values = {yes,no};
Interpolate.val = {yes};
Interpolate.help = {'Interpolate bad channels'};

S = cfg_exbranch;
S.tag = 'eeg_convertmat2nifti3D';
S.name = 'M/EEG convert2nifti3D';
S.val = {Fname n pixsize Interpolate};
S.help = {'Convert SPM M/EEG data to voxel-based images, as a time-series of 2D images'};
S.prog = @eeg_convertmat2nifti3D;
S.modality = {'EEG'};


function out = eeg_convertmat2nifti3D(job)
% construct the S struct
S = job;
S.Fname = strvcat(job.Fname);
if isfield(S.interpolate_bad, 'yes')
    S.interpolate_bad = 1;
else
    S.interpolate_bad = 0;
end

spm_eeg_convertmat2nifti3D(S);

