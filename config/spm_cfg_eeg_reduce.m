function S = spm_cfg_eeg_reduce
% Configuration file for M/EEG time-frequency analysis
%__________________________________________________________________________
% Copyright (C) 2010-2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_reduce.m 4911 2012-09-07 15:21:40Z vladimir $


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
% method
%--------------------------------------------------------------------------
method      = cfg_choice;
method.tag  = 'method';
method.name = 'Reduction method ';

specest_funs = spm_select('List',spm('dir'),'^spm_eeg_reduce_.*\.m$');
specest_funs = cellstr(specest_funs);
for i = 1:numel(specest_funs)
    method.values{i} = feval(spm_file(specest_funs{i},'basename'));
end

%--------------------------------------------------------------------------
% M/EEG Time-Frequency Analysis
%--------------------------------------------------------------------------
S = cfg_exbranch;
S.tag = 'analysis';
S.name = 'M/EEG Data reduction';
S.val = {D, spm_cfg_eeg_channel_selector, method};
S.help = {'Perform data reduction.'};
S.prog = @eeg_reduce;
S.vout = @vout_eeg_reduce;
S.modality = {'EEG'};

%==========================================================================
% function out = eeg_reduce(job)
%==========================================================================
function out = eeg_reduce(job)
% construct the S struct
S   = [];
S.D = job.D{1};

S.channels = spm_cfg_eeg_channel_selector(job.channels);

S.method = cell2mat(fieldnames(job.method));
S.settings = getfield(job.method, S.method);

D = spm_eeg_reduce(S);

out.D = D;
out.Dfname = {D.fname};


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