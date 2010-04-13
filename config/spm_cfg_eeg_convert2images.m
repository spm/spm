function S = spm_cfg_eeg_convert2images
% configuration file for writing voxel-based images from SPM M/EEG format,
% as a time-series of 2Dimages
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_convert2images.m 3818 2010-04-13 14:36:31Z vladimir $

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
S.tag = 'convert2images';
S.name = 'M/EEG Convert2Images';
S.val = {Fname n Interpolate};
S.help = {'Convert SPM M/EEG data to voxel-based images, as a time-series of 2D images'};
S.prog = @eeg_convert2images;
S.modality = {'EEG'};


function out = eeg_convert2images(job)
% construct the S struct
S = job;
S.Fname = strvcat(job.Fname);
if isfield(S.interpolate_bad, 'yes')
    S.interpolate_bad = 1;
else
    S.interpolate_bad = 0;
end

spm_eeg_convert2scalp(S);

