function reduce = spm_cfg_eeg_reduce
% Configuration file for M/EEG time-frequency analysis
%__________________________________________________________________________
% Copyright (C) 2010-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_reduce.m 6929 2016-11-14 13:07:31Z guillaume $


%--------------------------------------------------------------------------
% D
%--------------------------------------------------------------------------
D        = cfg_files;
D.tag    = 'D';
D.name   = 'File Name';
D.filter = 'mat';
D.num    = [1 Inf];
D.help   = {'Select the M/EEG mat file(s).'};


%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
method      = cfg_choice;
method.tag  = 'method';
method.name = 'Reduction method';
method.help = {'Reduction method'};

specest_funs = spm_select('List',spm('dir'),'^spm_eeg_reduce_.*\.m$');
specest_funs = cellstr(specest_funs);
method.values = cell(1,numel(specest_funs));
for i = 1:numel(specest_funs)
    method.values{i} = feval(spm_file(specest_funs{i},'basename'));
end

keeporig = cfg_menu;
keeporig.tag = 'keeporig';
keeporig.name = 'Keep original channels';
keeporig.labels = {'Yes', 'No'};
keeporig.values = {true, false};
keeporig.val = {false};
keeporig.help = {'Specify whether you want to keep the original unreduced channels'};

keepothers = cfg_menu;
keepothers.tag = 'keepothers';
keepothers.name = 'Keep other channels';
keepothers.labels = {'Yes', 'No'};
keepothers.values = {true, false};
keepothers.val = {false};
keepothers.help = {'Specify whether you want to keep channels that are not contributing to the new channels'};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''R''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'R'};

%--------------------------------------------------------------------------
% M/EEG Time-Frequency Analysis
%--------------------------------------------------------------------------
reduce = cfg_exbranch;
reduce.tag = 'reduce';
reduce.name = 'Data reduction';
reduce.val = {D, spm_cfg_eeg_channel_selector, method, keeporig, keepothers, prefix};
reduce.help = {'Perform data reduction.'};
reduce.prog = @eeg_reduce;
reduce.vout = @vout_eeg_reduce;
reduce.modality = {'EEG'};

%==========================================================================
% function out = eeg_reduce(job)
%==========================================================================
function out = eeg_reduce(job)
% construct the S struct
S   = [];
S.D = job.D{1};

S.prefix    = job.prefix;
S.channels  = spm_cfg_eeg_channel_selector(job.channels);

S.method    = cell2mat(fieldnames(job.method));
S.settings  = job.method.(S.method);
S.keepothers = job.keepothers;
S.keeporig  = job.keeporig;

D = spm_eeg_reduce(S);

out.D = D;
out.Dfname = {D.fullfile};


%==========================================================================
% function dep = vout_eeg_tf(job)
%==========================================================================
function dep = vout_eeg_reduce(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Reduced M/EEG dataset';
dep(1).src_output = substruct('.','D');
% this can be entered into any evaluated input
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Reduced M/EEG dataset';
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
