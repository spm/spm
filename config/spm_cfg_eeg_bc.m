function S = spm_cfg_eeg_bc
% configuration file for baseline correction
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_bc.m 3818 2010-04-13 14:36:31Z vladimir $

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
% time
%--------------------------------------------------------------------------
time         = cfg_entry;
time.tag     = 'time';
time.name    = 'Baseline';
time.help    = {'Start and stop of baseline [ms].'};
time.strtype = 'e';
time.num     = [1 2];

%--------------------------------------------------------------------------
% S
%--------------------------------------------------------------------------
S          = cfg_exbranch;
S.tag      = 'bc';
S.name     = 'M/EEG Baseline correction';
S.val      = {D, time};
S.help     = {'Baseline correction of M/EEG time data'}';
S.prog     = @eeg_bc;
S.vout     = @vout_eeg_bc;
S.modality = {'EEG'};

%==========================================================================
function out = eeg_bc(job)
% construct the S struct
S.D            = job.D{1};
S.time         = job.time;

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
