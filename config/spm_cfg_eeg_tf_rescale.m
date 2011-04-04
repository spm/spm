function S = spm_cfg_eeg_tf_rescale
% configuration file for rescaling spectrograms
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_cfg_eeg_tf_rescale.m 4287 2011-04-04 13:55:54Z vladimir $

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
% Db
%--------------------------------------------------------------------------
Db        = cfg_files;
Db.tag    = 'Db';
Db.name   = 'External baseline dataset';
Db.filter = 'mat';
Db.num    = [0 1];
Db.val    = {[]};
Db.help   = {'Select the baseline M/EEG mat file. Leave empty to use the input dataset'};

%--------------------------------------------------------------------------
% Sbaseline
%--------------------------------------------------------------------------
Sbaseline         = cfg_entry;
Sbaseline.tag     = 'Sbaseline';
Sbaseline.name    = 'Baseline';
Sbaseline.help    = {'Start and stop of baseline [ms].'};
Sbaseline.strtype = 'e';
Sbaseline.num     = [1 2];

%--------------------------------------------------------------------------
% baseline
%--------------------------------------------------------------------------
baseline         = cfg_branch;
baseline.tag     = 'baseline';
baseline.name    = 'Baseline';
baseline.help    = {'Baseline parameters.'};
baseline.val     = {Sbaseline, Db};

%--------------------------------------------------------------------------
% method_logr
%--------------------------------------------------------------------------
method_logr      = cfg_branch;
method_logr.tag  = 'LogR';
method_logr.name = 'Log Ratio';
method_logr.val  = {baseline};
method_logr.help = {'Log Ratio.'};

%--------------------------------------------------------------------------
% method_diff
%--------------------------------------------------------------------------
method_diff      = cfg_branch;
method_diff.tag  = 'Diff';
method_diff.name = 'Difference';
method_diff.val  = {baseline};
method_diff.help = {'Difference.'};

%--------------------------------------------------------------------------
% method_rel
%--------------------------------------------------------------------------
method_rel       = cfg_branch;
method_rel.tag   = 'Rel';
method_rel.name  = 'Relative';
method_rel.val   = {baseline};
method_rel.help  = {'Relative.'};

%--------------------------------------------------------------------------
% method_zscore
%--------------------------------------------------------------------------
method_zscore       = cfg_branch;
method_zscore.tag   = 'Zscore';
method_zscore.name  = 'Zscore';
method_zscore.val   = {baseline};
method_zscore.help  = {'Z score'};

%--------------------------------------------------------------------------
% method_log
%--------------------------------------------------------------------------
method_log       = cfg_const;
method_log.tag   = 'Log';
method_log.name  = 'Log';
method_log.val   = {1};
method_log.help  = {'Log.'};

%--------------------------------------------------------------------------
% method_sqrt
%--------------------------------------------------------------------------
method_sqrt      = cfg_const;
method_sqrt.tag  = 'Sqrt';
method_sqrt.name = 'Sqrt';
method_sqrt.val  = {1};
method_sqrt.help = {'Square Root.'};

%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
method        = cfg_choice;
method.tag    = 'method';
method.name   = 'Rescale method';
method.val    = {method_logr};
method.help   = {'Select the rescale method.'};
method.values = {method_logr method_diff method_rel method_zscore method_log method_sqrt};

%--------------------------------------------------------------------------
% S
%--------------------------------------------------------------------------
S          = cfg_exbranch;
S.tag      = 'rescale';
S.name     = 'M/EEG Time-Frequency Rescale';
S.val      = {D, method};
S.help     = {'Rescale (avg) spectrogram with nonlinear and/or difference operator.'
              'For ''Log'' and ''Sqrt'', these functions are applied to spectrogram.'
              'For ''LogR'', ''Rel'' and ''Diff'' this function computes power in the baseline.'
              'p_b and outputs:'
              '(i) p-p_b for ''Diff'''
              '(ii) 100*(p-p_b)/p_b for ''Rel'''
              '(iii) log (p/p_b) for ''LogR'''}';
S.prog     = @eeg_tf_rescale;
S.vout     = @vout_eeg_tf_rescale;
S.modality = {'EEG'};

%==========================================================================
function out = eeg_tf_rescale(job)
% construct the S struct
S.D            = job.D{1};
S.tf.method    = fieldnames(job.method);
S.tf.method    = S.tf.method{1};
switch lower(S.tf.method)
    case {'logr','diff', 'rel', 'zscore'}
        S.tf.Sbaseline = 1e-3*job.method.(S.tf.method).baseline.Sbaseline;
        if ~(isempty(job.method.(S.tf.method).baseline.Db) || isequal(job.method.(S.tf.method).baseline.Db, {''}))
            S.tf.Db = job.method.(S.tf.method).baseline.Db{1};
        end
    case {'log', 'sqrt'}
end

out.D          = spm_eeg_tf_rescale(S);
out.Dfname     = {fullfile(out.D.path,out.D.fname)};

%==========================================================================
function dep = vout_eeg_tf_rescale(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Rescaled TF Data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Rescaled TF Datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
