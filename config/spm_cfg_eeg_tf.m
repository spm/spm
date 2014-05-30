function tf = spm_cfg_eeg_tf
% Configuration file for M/EEG time-frequency analysis
%__________________________________________________________________________
% Copyright (C) 2010-2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_tf.m 6029 2014-05-30 18:52:03Z vladimir $


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
timewin.strtype = 'r';
timewin.num     = [1 2];
timewin.val     = {[-Inf Inf]};
timewin.help    = {'Time window (ms)'};

%--------------------------------------------------------------------------
% frequencies
%--------------------------------------------------------------------------
frequencies         = cfg_entry;
frequencies.tag     = 'frequencies';
frequencies.name    = 'Frequencies of interest';
frequencies.strtype = 'r';
frequencies.num     = [0 Inf];
frequencies.val     = {[]};
frequencies.help    = {'Frequencies of interest (as a vector), if empty 1-48 with optimal frequency bins ~1 Hz or above resolution'};

%--------------------------------------------------------------------------
% phase
%--------------------------------------------------------------------------
phase        = cfg_menu;
phase.tag    = 'phase';
phase.name   = 'Save phase';
phase.help   = {'Save phase as well as power'};
phase.labels = {'yes', 'no'};
phase.values = {1, 0};
phase.val    = {0};

%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
method      = cfg_choice;
method.tag  = 'method';
method.name = 'Spectral estimation ';

specest_funs = spm_select('List',spm('dir'),'^spm_eeg_specest_.*\.m$');
specest_funs = cellstr(specest_funs);
for i = 1:numel(specest_funs)
    method.values{i} = feval(spm_file(specest_funs{i},'basename'));
end

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset.'};
prefix.strtype = 's';
prefix.num     = [0 Inf];
prefix.val     = {''};

%--------------------------------------------------------------------------
% M/EEG Time-Frequency Analysis
%--------------------------------------------------------------------------
tf = cfg_exbranch;
tf.tag = 'tf';
tf.name = 'Time-frequency analysis';
tf.val = {D, spm_cfg_eeg_channel_selector, frequencies, timewin, method, phase, prefix};
tf.help = {'Perform time-frequency analysis of epoched M/EEG data.'};
tf.prog = @eeg_tf;
tf.vout = @vout_eeg_tf;
tf.modality = {'EEG'};

%==========================================================================
% function out = eeg_tf(job)
%==========================================================================
function out = eeg_tf(job)
% construct the S struct
S   = [];
S.D = job.D{1};

S.channels = spm_cfg_eeg_channel_selector(job.channels);

S.frequencies = job.frequencies;
S.timewin = job.timewin;
S.phase = job.phase;

S.method = cell2mat(fieldnames(job.method));
S.settings = job.method.(S.method);

S.prefix = job.prefix;

[Dtf, Dtph] = spm_eeg_tf(S);

out.Dtf = Dtf;
out.Dtfname = {Dtf.fullfile};

out.Dtph = Dtph;
if ~isempty(Dtph)
    out.Dtphname = {Dtph.fullfile};
else
    out.Dtphname = {''};
end

%==========================================================================
% function dep = vout_eeg_tf(job)
%==========================================================================
function dep = vout_eeg_tf(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'M/EEG time-frequency power dataset';
dep(1).src_output = substruct('.','Dtf');
% this can be entered into any evaluated input
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'M/EEG time-frequency power dataset';
dep(2).src_output = substruct('.','Dtfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});

if job.phase
    dep(3)            = cfg_dep;
    dep(3).sname      = 'M/EEG time-frequency phase dataset';
    dep(3).src_output = substruct('.','Dtph');
    % this can be entered into any evaluated input
    dep(3).tgt_spec   = cfg_findspec({{'strtype','e'}});
    
    dep(4)            = cfg_dep;
    dep(4).sname      = 'M/EEG time-frequency phase dataset';
    dep(4).src_output = substruct('.','Dtphname');
    % this can be entered into any file selector
    dep(4).tgt_spec   = cfg_findspec({{'filter','mat'}});
end