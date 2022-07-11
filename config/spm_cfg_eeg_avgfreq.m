function avgfreq = spm_cfg_eeg_avgfreq
% Configuration file for averaging over frequency
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%--------------------------------------------------------------------------
% avgfreq
%--------------------------------------------------------------------------
avgfreq          = cfg_exbranch;
avgfreq.tag      = 'avgfreq';
avgfreq.name     = 'Average over frequency';
avgfreq.val      = @avgfreq_cfg;
avgfreq.help     = {'Average M/EEG data over frequency'}';
avgfreq.prog     = @eeg_avgfreq;
avgfreq.vout     = @vout_eeg_avgfreq;
avgfreq.modality = {'EEG'};


%==========================================================================
function varargout = avgfreq_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

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
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''P''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'P'};

[cfg,varargout{1}] = deal({D, freqwin, prefix});


%==========================================================================
function out = eeg_avgfreq(job)
% construct the S struct
S           = job;
S.D         = S.D{1};
out.D       = spm_eeg_avgfreq(S);
out.Dfname  = {fullfile(out.D)};

%==========================================================================
function dep = vout_eeg_avgfreq(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Frequency averaged M/EEG data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Frequency averaged M/EEG datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
