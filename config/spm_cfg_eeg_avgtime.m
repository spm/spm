function avgtime = spm_cfg_eeg_avgtime
% configuration file for averaging over time
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_avgtime.m 5652 2013-09-25 09:36:22Z volkmar $

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
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''S''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'S'};

%--------------------------------------------------------------------------
% avgfreq
%--------------------------------------------------------------------------
avgtime          = cfg_exbranch;
avgtime.tag      = 'avgtime';
avgtime.name     = 'Average over time';
avgtime.val      = {D, timewin, prefix};
avgtime.help     = {'Average M/EEG data over time'}';
avgtime.prog     = @eeg_avgtime;
avgtime.vout     = @vout_eeg_avgtime;
avgtime.modality = {'EEG'};

%==========================================================================
function out = eeg_avgtime(job)
% construct the S struct
S           = job;
S.D         = S.D{1};
out.D       = spm_eeg_avgtime(S);
out.Dfname  = {fullfile(out.D)};

%==========================================================================
function dep = vout_eeg_avgtime(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Time averaged M/EEG data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Time averaged M/EEG datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
