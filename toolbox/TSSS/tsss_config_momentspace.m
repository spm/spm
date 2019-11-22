function momentspace = tsss_config_momentspace
% configuration file for cropping
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: tsss_config_momentspace.m 7703 2019-11-22 12:06:29Z guillaume $

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
% condthresh
%--------------------------------------------------------------------------
condthresh         = cfg_entry;
condthresh.tag     = 'condthresh';
condthresh.name    = 'Condition number threshold';
condthresh.help    = {'Threshold on condition number applied for basis regularisation'};
condthresh.strtype = 'r';
condthresh.num     = [1 1];
condthresh.val     = {80};

%--------------------------------------------------------------------------
% none
%--------------------------------------------------------------------------

none = cfg_const;
none.tag = 'none';
none.name = 'None';
none.val  = {0};

%--------------------------------------------------------------------------
% addchannels
%--------------------------------------------------------------------------

addchannels      = cfg_choice;
addchannels.tag  = 'addchannels';
addchannels.name = 'Extra channels to add';
addchannels.values  = {none, spm_cfg_eeg_channel_selector};
addchannels.val  = {none};

%--------------------------------------------------------------------------
% momentspace
%--------------------------------------------------------------------------
momentspace          = cfg_exbranch;
momentspace.tag      = 'momentspace';
momentspace.name     = 'TSSS space conversion';
momentspace.val      = {D, condthresh, addchannels};
momentspace.help     = {'Engage the TSSS space virtual montage'}';
momentspace.prog     = @eeg_momentspace;
momentspace.vout     = @vout_eeg_momentspace;
momentspace.modality = {'EEG'};

%==========================================================================
function out = eeg_momentspace(job)
% construct the S struct
S   = [];
S.D = char(job.D);
S.condthresh = job.condthresh;

D = spm_eeg_load(S.D);

if isfield(job.addchannels, 'channels')
    S.addchannels = D.chanlabels(D.selectchannels(spm_cfg_eeg_channel_selector(job.addchannels.channels)));
else
    S.addchannels = {};
end

out.D        = tsss_spm_momentspace(S);
out.Dfname   = {fullfile(out.D)};

%==========================================================================
function dep = vout_eeg_momentspace(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Moment space MEG data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Moment space MEG datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});