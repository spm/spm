function S = spm_cfg_eeg_grandmean
% configuration file for averaging evoked responses
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_grandmean.m 3881 2010-05-07 21:02:57Z vladimir $

% -------------------------------------------------------------------------
% weighted Weighted average
% -------------------------------------------------------------------------
weighted         = cfg_menu;
weighted.tag     = 'weighted';
weighted.name    = 'Weighted average?';
weighted.help    = {'Average weighted by number of replications in input.'};
weighted.labels  = {'Yes' 'No'};
weighted.values  = {1 0};

%--------------------------------------------------------------------------
% D File Names
%--------------------------------------------------------------------------
D            = cfg_files;
D.tag        = 'D';
D.name       = 'File Names';
D.filter     = 'mat';
D.num        = [1 inf];
D.help       = {'Select the M/EEG mat file.'};

%--------------------------------------------------------------------------
% Dout Output filename
%--------------------------------------------------------------------------
Dout         = cfg_entry;
Dout.tag     = 'Dout';
Dout.name    = 'Output filename';
Dout.strtype = 's';
Dout.num     = [1 inf];
Dout.help    = {'Choose filename'};

%--------------------------------------------------------------------------
% S Grandmean
%--------------------------------------------------------------------------
S            = cfg_exbranch;
S.tag        = 'grandmean';
S.name       = 'M/EEG Grandmean';
S.val        = {D Dout weighted};
S.help       = {'Average multiple evoked responses'};
S.prog       = @eeg_grandmean;
S.vout       = @vout_eeg_grandmean;
S.modality   = {'EEG'};

%==========================================================================
function out = eeg_grandmean(job)
% construct the S struct
S.D = strvcat(job.D);
S.Dout = job.Dout;
S.weighted = job.weighted;

out.D = spm_eeg_grandmean(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

%==========================================================================
function dep = vout_eeg_grandmean(job)

dep(1)            = cfg_dep;
dep(1).sname      = 'Grandmean Data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Grandmean Datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
