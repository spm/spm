function bc = spm_cfg_eeg_bc
% configuration file for baseline correction
%__________________________________________________________________________
% Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_bc.m 5652 2013-09-25 09:36:22Z volkmar $

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
timewin.name    = 'Baseline';
timewin.help    = {'Start and stop of baseline [ms].'};
timewin.strtype = 'r';
timewin.num     = [1 2];

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the filtered dataset. Default prefix is ''b''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'b'};

%--------------------------------------------------------------------------
% bc
%--------------------------------------------------------------------------
bc          = cfg_exbranch;
bc.tag      = 'bc';
bc.name     = 'Baseline correction';
bc.val      = {D, timewin, prefix};
bc.help     = {'Baseline correction of M/EEG time data'}';
bc.prog     = @eeg_bc;
bc.vout     = @vout_eeg_bc;
bc.modality = {'EEG'};

%==========================================================================
function out = eeg_bc(job)
% construct the S struct
S = job;
S.D = S.D{1};

out.D          = spm_eeg_bc(S);
out.Dfname     = {fullfile(out.D.path,out.D.fname)};

%==========================================================================
function dep = vout_eeg_bc(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Baseline corrected M/EEG data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Baseline corrected M/EEG datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
