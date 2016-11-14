function cfc = spm_cfg_eeg_cfc
% Configuration file for M/EEG cross-frequency coupling analysis
%__________________________________________________________________________
% Copyright (C) 2014-2016 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_cfc.m 6929 2016-11-14 13:07:31Z guillaume $


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
% window length
%--------------------------------------------------------------------------
window        = cfg_entry;
window.tag     = 'window';
window.name    = 'Window length';
window.help    = {'Time window length (ms). Only used for continuous data.'};
window.strtype = 'r';
window.num     = [1 1];
window.val     = {1000};

%--------------------------------------------------------------------------
% regressors of interest
%--------------------------------------------------------------------------
regressors      = cfg_repeat;
regressors.tag  = 'regressors';
regressors.name = 'Regressors of interest';
regressors.num  = [1 Inf];
regressors.help = {'Regressors of interest'};

reg_funs = {'spm_eeg_regressors_tfpower.m', 'spm_eeg_regressors_tfphase.m'};
regressors.values = cell(1,numel(reg_funs));
for i = 1:numel(reg_funs)
    regressors.values{i} = feval(spm_file(reg_funs{i},'basename'));
end

%--------------------------------------------------------------------------
% confounds
%--------------------------------------------------------------------------
confounds      = cfg_repeat;
confounds.tag  = 'confounds';
confounds.name = 'Confounds';
confounds.num  = [0 Inf];
confounds.help = {'Confounds'};

reg_funs = spm_select('List',spm('dir'),'^spm_eeg_regressors_.*\.m$');
reg_funs = cellstr(reg_funs);
confounds.values = cell(1,numel(reg_funs));
for i = 1:numel(reg_funs)
    confounds.values{i} = feval(spm_file(reg_funs{i},'basename'));
end

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Directory prefix';
prefix.help    = {'Specify the string to be prepended to the output directory name'};
prefix.strtype = 's';
prefix.num     = [0 Inf];
prefix.val     = {''};

%--------------------------------------------------------------------------
% M/EEG Cross-Frequency Coupling Analysis
%--------------------------------------------------------------------------
cfc = cfg_exbranch;
cfc.tag = 'cfc';
cfc.name = 'Cross-frequency coupling';
cfc.val = {D, spm_cfg_eeg_channel_selector, conditions, freqwin, window, regressors, confounds, prefix};
cfc.help = {'GLM-based cross-frequency coupling analysis'};
cfc.prog = @eeg_cfc;
%cfc.vout = @vout_eeg_cfc;
cfc.modality = {'EEG'};

%==========================================================================
% function out = eeg_cfc(job)
%==========================================================================
function out = eeg_cfc(job)
S          = job;
S.D        = char(job.D);
S.channels = spm_cfg_eeg_channel_selector(job.channels);
spm_eeg_cfc(S);
out = [];
% out.regrfile = {spm_eeg_regressors(job)};
% out.inputfile = {job.D};

%==========================================================================
% function dep = vout_eeg_tf(job)
%==========================================================================
% function dep = vout_eeg_cfc(job)
% % return dependencies
% dep(1)            = cfg_dep;
% dep(1).sname      = 'M/EEG dataset for statistics';
% dep(1).src_output = substruct('.','inputfile');
% dep(1).tgt_spec   = cfg_findspec({{'filter','mat'}});
% dep(2)            = cfg_dep;
% dep(2).sname      = 'MEEG GLM regressors';
% dep(2).src_output = substruct('.','regrfile');
% dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
% 
% 
