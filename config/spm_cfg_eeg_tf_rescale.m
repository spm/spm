function rescale = spm_cfg_eeg_tf_rescale
% configuration file for rescaling spectrograms
%__________________________________________________________________________
% Copyright (C) 2009-2013 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_cfg_eeg_tf_rescale.m 6825 2016-07-04 10:03:57Z vladimir $

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
% timewin
%--------------------------------------------------------------------------
timewin         = cfg_entry;
timewin.tag     = 'timewin';
timewin.name    = 'Baseline time window';
timewin.help    = {'Start and stop of the baseline time window [ms].'};
timewin.strtype = 'r';
timewin.num     = [1 2];
timewin.val     = {[-Inf 0]};

%--------------------------------------------------------------------------
% pooledbaseline
%--------------------------------------------------------------------------

pooledbaseline = cfg_menu;
pooledbaseline.tag = 'pooledbaseline';
pooledbaseline.name = 'Pool baseline across trials';
pooledbaseline.labels = {'Yes', 'No'};
pooledbaseline.val = {0};
pooledbaseline.values = {1,0};
pooledbaseline.help = {'Combine baseline across trials to avoid bias',...
    'See Ciuparu and Muresan Eur J Neurosci. 43(7):861-9, 2016'};

%--------------------------------------------------------------------------
% baseline
%--------------------------------------------------------------------------
baseline         = cfg_branch;
baseline.tag     = 'baseline';
baseline.name    = 'Baseline';
baseline.help    = {'Baseline parameters.'};
baseline.val     = {timewin, pooledbaseline, Db};

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
% method_logeps
%--------------------------------------------------------------------------
method_logeps       = cfg_const;
method_logeps.tag   = 'Logeps';
method_logeps.name  = 'Log+eps';
method_logeps.val   = {1};
method_logeps.help  = {'Log + epsilon (to avoid boosting very low values)'};

%--------------------------------------------------------------------------
% method_sqrt
%--------------------------------------------------------------------------
method_sqrt      = cfg_const;
method_sqrt.tag  = 'Sqrt';
method_sqrt.name = 'Sqrt';
method_sqrt.val  = {1};
method_sqrt.help = {'Square Root.'};


%--------------------------------------------------------------------------
% method_sqrt
%--------------------------------------------------------------------------
method_none      = cfg_const;
method_none.tag  = 'None';
method_none.name = 'None';
method_none.val  = {1};
method_none.help = {'No rescaling - just copy. Useful for optimising batch pipelines.'};

%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
method        = cfg_choice;
method.tag    = 'method';
method.name   = 'Rescale method';
method.val    = {method_logr};
method.help   = {'Select the rescale method.'};
method.values = {method_logr method_diff method_rel method_zscore method_log method_logeps method_sqrt, method_none};

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''r''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'r'};

%--------------------------------------------------------------------------
% rescale
%--------------------------------------------------------------------------
rescale          = cfg_exbranch;
rescale.tag      = 'rescale';
rescale.name     = 'Time-frequency rescale';
rescale.val      = {D, method, prefix};
rescale.help     = {'Rescale (avg) spectrogram with nonlinear and/or difference operator.'
              'For ''Log'' and ''Sqrt'', these functions are applied to spectrogram.'
              'For ''LogR'', ''Rel'' and ''Diff'' this function computes power in the baseline.'
              'p_b and outputs:'
              '(i) p-p_b for ''Diff'''
              '(ii) 100*(p-p_b)/p_b for ''Rel'''
              '(iii) log (p/p_b) for ''LogR'''}';
rescale.prog     = @eeg_tf_rescale;
rescale.vout     = @vout_eeg_tf_rescale;
rescale.modality = {'EEG'};

%==========================================================================
function out = eeg_tf_rescale(job)
% construct the S struct
S.D         = job.D{1};
S.method    = char(fieldnames(job.method));
S.prefix    = job.prefix;

if ismember(lower(S.method), {'logr','diff', 'rel', 'zscore'})
    S.timewin        = job.method.(S.method).baseline.timewin;
    S.pooledbaseline = job.method.(S.method).baseline.pooledbaseline;
    if ~(isempty(job.method.(S.method).baseline.Db) || isequal(job.method.(S.method).baseline.Db, {''}))
        S.Db = job.method.(S.method).baseline.Db{1};
    end
end

out.D          = spm_eeg_tf_rescale(S);
out.Dfname     = {fullfile(out.D)};

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
