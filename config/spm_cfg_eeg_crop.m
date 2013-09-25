function crop = spm_cfg_eeg_crop
% configuration file for cropping
%__________________________________________________________________________
% Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_crop.m 5652 2013-09-25 09:36:22Z volkmar $

%--------------------------------------------------------------------------
% D
%--------------------------------------------------------------------------
D        = cfg_files;
D.tag    = 'D';
D.name   = 'File Name';
D.filter = 'mat';
D.num    = [1 1];
D.help   = {'Select the M/EEG mat file.'};

%--------------------------------------------------------------------------
% timewin
%--------------------------------------------------------------------------
timewin         = cfg_entry;
timewin.tag     = 'timewin';
timewin.name    = 'Time window';
timewin.help    = {'Start and stop of the time window [ms].'};
timewin.strtype = 'r';
timewin.num     = [1 2];
timewin.val     = {[-Inf Inf]};

%--------------------------------------------------------------------------
% freqwin
%--------------------------------------------------------------------------
freqwin         = cfg_entry;
freqwin.tag     = 'freqwin';
freqwin.name    = 'Frequency window';
freqwin.help    = {'Start and stop of the frequency window (Hz).'};
freqwin.strtype = 'r';
freqwin.num     = [1 2];
freqwin.val     = {[-Inf Inf]};

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''p''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'p'};

%--------------------------------------------------------------------------
% crop
%--------------------------------------------------------------------------
crop          = cfg_exbranch;
crop.tag      = 'crop';
crop.name     = 'Crop';
crop.val      = {D, timewin, freqwin, spm_cfg_eeg_channel_selector, prefix};
crop.help     = {'Cropping M/EEG data'}';
crop.prog     = @eeg_crop;
crop.vout     = @vout_eeg_crop;
crop.modality = {'EEG'};

%==========================================================================
function out = eeg_crop(job)
% construct the S struct
S           = job;
S.D         = S.D{1};
S.channels  = spm_cfg_eeg_channel_selector(job.channels);
out.D       = spm_eeg_crop(S);
out.Dfname  = {fullfile(out.D)};

%==========================================================================
function dep = vout_eeg_crop(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Cropped M/EEG data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Cropped M/EEG datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
