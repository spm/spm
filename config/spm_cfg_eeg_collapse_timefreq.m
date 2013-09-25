function collapse = spm_cfg_eeg_collapse_timefreq
% configuration file for within-image averaging
%__________________________________________________________________________
% Copyright (C) 2009-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_collapse_timefreq.m 5652 2013-09-25 09:36:22Z volkmar $

% ---------------------------------------------------------------------
% Images to Smooth
% ---------------------------------------------------------------------
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Images to average';
images.help    = {'Specify the images to average time/frequency.'};
images.filter = 'image';
images.ufilter = '.*';
images.num     = [0 Inf];

%--------------------------------------------------------------------------
% timewin
%--------------------------------------------------------------------------
timewin         = cfg_entry;
timewin.tag     = 'timewin';
timewin.name    = 'Time/frequency window';
timewin.help    = {'Start and stop of the time/frequency window [ms/Hz].'};
timewin.strtype = 'r';
timewin.num     = [1 2];
timewin.val     = {[-Inf Inf]};

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the output images. Default prefix is ''l''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'l'};

%--------------------------------------------------------------------------
% crop
%--------------------------------------------------------------------------
collapse          = cfg_exbranch;
collapse.tag      = 'collapse';
collapse.name     = 'Collapse time';
collapse.val      = {images, timewin, prefix};
collapse.help     = {'Compute within-peristimulus time (or frequency) averages (contrasts) of M/EEG data in voxel-space'}';
collapse.prog     = @eeg_collapse;
collapse.vout     = @vout_eeg_collapse;
collapse.modality = {'EEG'};

%==========================================================================
function out = eeg_collapse(job)
out.files   = spm_eeg_collapse_timefreq(job);

%==========================================================================
function dep = vout_eeg_collapse(job)
% Output file names will be saved in a struct with field .files
dep(1)            = cfg_dep;
dep(1).sname      = 'Collapsed M/EEG images';
dep(1).src_output = substruct('.','files');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

