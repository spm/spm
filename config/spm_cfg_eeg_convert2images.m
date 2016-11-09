function convert2images = spm_cfg_eeg_convert2images
% Configuration file for writing voxel-based images from SPM M/EEG format,
% as a time-series of 2Dimages
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_convert2images.m 6926 2016-11-09 22:13:19Z guillaume $

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
% mode
%--------------------------------------------------------------------------
mode = cfg_menu;
mode.tag = 'mode';
mode.name = 'Mode';
mode.labels = {
    'scalp x time'
    'scalp x frequency' 
    'scalp' 
    'source'
    'time x frequency' 
    'time' 
    'frequency'
    'average'}';
mode.values = mode.labels;
mode.help = {'Select the mode for conversion to images.'};

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
% conditions
%--------------------------------------------------------------------------
condlabel = cfg_entry;
condlabel.tag = 'conditions';
condlabel.name = 'Condition label';
condlabel.strtype = 's';
condlabel.help = {''};

conditions = cfg_repeat;
conditions.tag = 'condrepeat';
conditions.name = 'Conditions';
conditions.help = {'Specify the labels of the conditions to be converted.'};
conditions.num  = [0 Inf];
conditions.values  = {condlabel};
conditions.val = {};

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Directory prefix';
prefix.help    = {'Specify the string to be prepended to the output directory name.'};
prefix.strtype = 's';
prefix.num     = [0 Inf];
prefix.val     = {''};


convert2images = cfg_exbranch;
convert2images.tag = 'convert2images';
convert2images.name = 'Convert2Images';
convert2images.val = {D, mode, conditions, spm_cfg_eeg_channel_selector, timewin, freqwin, prefix};
convert2images.help = {'Convert SPM M/EEG data to voxel-based images.'};
convert2images.prog = @run_convert2images;
convert2images.vout = @vout_convert2images;
convert2images.modality = {'EEG'};


%==========================================================================
function out = run_convert2images(job)

S           = job;
S.D         = S.D{1};
S.channels  = spm_cfg_eeg_channel_selector(job.channels);

[out.files, out.dir{1}] = spm_eeg_convert2images(S);


%==========================================================================
function dep = vout_convert2images(varargin)
% Output file names will be saved in a struct with field .files
dep(1)            = cfg_dep;
dep(1).sname      = 'M/EEG exported images';
dep(1).src_output = substruct('.','files');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'M/EEG images folder';
dep(2).src_output = substruct('.','dir');
dep(2).tgt_spec   = cfg_findspec({{'filter','dir','strtype','e'}});
