function copy = spm_cfg_eeg_copy
% Configuration file for copying M/EEG datasets
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2009-2022 Wellcome Centre for Human Neuroimaging


%--------------------------------------------------------------------------
% copy
%--------------------------------------------------------------------------
copy          = cfg_exbranch;
copy.tag      = 'copy';
copy.name     = 'Copy';
copy.val      = @copy_cfg;
copy.help     = {'Copying M/EEG datasets'}';
copy.prog     = @eeg_copy;
copy.vout     = @vout_eeg_copy;
copy.modality = {'EEG'};


%==========================================================================
function varargout = copy_cfg

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
% outfile
%--------------------------------------------------------------------------
outfile         = cfg_entry;
outfile.tag     = 'outfile';
outfile.name    = 'Output filename';
outfile.strtype = 's';
outfile.num     = [0 inf];
outfile.help    = {'Choose filename.'};

[cfg,varargout{1}] = deal({D, outfile});


%==========================================================================
function out = eeg_copy(job)
% construct the S struct
S           = job;
S.D         = S.D{1};
out.D       = spm_eeg_copy(S);
out.Dfname  = {fullfile(out.D)};

%==========================================================================
function dep = vout_eeg_copy(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Copied M/EEG data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Copied M/EEG datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
