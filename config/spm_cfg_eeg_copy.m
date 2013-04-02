function copy = spm_cfg_eeg_copy
% configuration file for copying
%__________________________________________________________________________
% Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_copy.m 5377 2013-04-02 17:07:57Z vladimir $

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
outfile = cfg_entry;
outfile.tag = 'outfile';
outfile.name = 'Output filename';
outfile.strtype = 's';
outfile.num = [0 inf];
outfile.help = {'Choose filename.'};

%--------------------------------------------------------------------------
% copy
%--------------------------------------------------------------------------
copy          = cfg_exbranch;
copy.tag      = 'copy';
copy.name     = 'Copy';
copy.val      = {D, outfile};
copy.help     = {'Copying M/EEG datasets'}';
copy.prog     = @eeg_copy;
copy.vout     = @vout_eeg_copy;
copy.modality = {'EEG'};

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
