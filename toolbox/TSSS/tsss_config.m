function tsss = tsss_config
% configuration file for cropping
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: tsss_config.m 7703 2019-11-22 12:06:29Z guillaume $

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
% temporal
%--------------------------------------------------------------------------

temporal = cfg_menu;
temporal.tag = 'temporal';
temporal.name = 'Do temporal denoising';
temporal.labels = {'Yes', 'No'};
temporal.val = {1};
temporal.values = {1,0};
temporal.help = {'Determines whether to use temporal denoising (tSSS)',...
    'if not, just SSS is done.'};


%--------------------------------------------------------------------------
% realign
%--------------------------------------------------------------------------
Dref        = cfg_files;
Dref.tag    = 'Dref';
Dref.name   = 'Reference dataset';
Dref.filter = 'mat';
Dref.num    = [0 1];
Dref.val    = {[]};
Dref.help   = {'Select the M/EEG mat file.',...
    'Leave empty to realign within dataset'};


refind         = cfg_entry;
refind.tag     = 'refind';
refind.name    = 'Reference index';
refind.help    = {'Index of the reference sensors within dataset',...
    'Normally should be 1. Set to 0 for not realigning a composite dataset'};
refind.strtype = 'w';
refind.num     = [1 1];
refind.val     = {1};


realign         = cfg_branch;
realign.tag      = 'realign';
realign.name     = 'Realign head location';
realign.val      = {Dref, refind};
realign.help     = {'Realign head location to another dataset'}';

%--------------------------------------------------------------------------
% timewin
%--------------------------------------------------------------------------
timewin         = cfg_entry;
timewin.tag     = 'timewin';
timewin.name    = 'Time window';
timewin.help    = {'Time window (in sec) for temporal correlation',...
    'Ignored for epoched data'};
timewin.strtype = 'r';
timewin.num     = [1 1];
timewin.val     = {1};

%--------------------------------------------------------------------------
% corrlimit
%--------------------------------------------------------------------------
corrlimit         = cfg_entry;
corrlimit.tag     = 'corrlimit';
corrlimit.name    = 'Correlation limit';
corrlimit.help    = {'Correlation limit parameter for tSSS'};
corrlimit.strtype = 'r';
corrlimit.num     = [1 1];
corrlimit.val     = {0.98};

%--------------------------------------------------------------------------
% Lin
%--------------------------------------------------------------------------
Lin         = cfg_entry;
Lin.tag     = 'Lin';
Lin.name    = 'Inner dimension';
Lin.help    = {'Order if the inner SSS basis'};
Lin.strtype = 'n';
Lin.num     = [1 1];
Lin.val     = {8};

%--------------------------------------------------------------------------
% Lout
%--------------------------------------------------------------------------
Lout         = cfg_entry;
Lout.tag     = 'Lout';
Lout.name    = 'Outer dimension';
Lout.help    = {'Order if the outer SSS basis'};
Lout.strtype = 'n';
Lout.num     = [1 1];
Lout.val     = {3};

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
% ospace
%--------------------------------------------------------------------------

ospace = cfg_menu;
ospace.tag = 'ospace';
ospace.name = 'Output space';
ospace.labels = {'sensor', 'SSS'};
ospace.val = {0};
ospace.values = {0,1};
ospace.help = {'Determines whether the output file is in sensor space',...
    'or has virtual montage trasforming to SSS space.'};

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''tsss_''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'tsss_'};

%--------------------------------------------------------------------------
% tsss
%--------------------------------------------------------------------------
tsss          = cfg_exbranch;
tsss.tag      = 'tsss';
tsss.name     = 'TSSS denoising';
tsss.val      = {D, temporal, realign, timewin, corrlimit, Lin, Lout, condthresh, ospace, prefix};
tsss.help     = {'TSSS clean-up for Neuromag data'}';
tsss.prog     = @eeg_tsss;
tsss.vout     = @vout_eeg_tsss;
tsss.modality = {'EEG'};

%==========================================================================
function out = eeg_tsss(job)
% construct the S struct
S = [];
S.D = char(job.D);
S.tsss       = job.temporal;
S.Dref       = char(job.realign.Dref);
S.refind     = job.realign.refind;
S.t_window   = job.timewin;
S.corr_limit = job.corrlimit;
S.Lin        = job.Lin;
S.Lout       = job.Lout;
S.cond_threshold = job.condthresh;
S.xspace     = job.ospace;
S.prefix     = job.prefix;
out.D        = tsss_spm_enm(S);
out.Dfname   = {fullfile(out.D)};

%==========================================================================
function dep = vout_eeg_tsss(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'TSSS-ed MEG data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'TSSS-ed MEG datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});