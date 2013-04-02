function combineplanar = spm_cfg_eeg_combineplanar
% configuration file for combineplanar
%__________________________________________________________________________
% Copyright (C) 2009-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_combineplanar.m 5377 2013-04-02 17:07:57Z vladimir $

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
% mode
%--------------------------------------------------------------------------
mode = cfg_menu;
mode.tag = 'mode';
mode.name = 'Copying mode';
mode.labels = {'Append', 'Replace planar', 'Replace all MEG', 'Only keep combined'};
mode.val = {'replace'};
mode.values = {'append', 'replace', 'replacemeg', 'keep'};
mode.help = {'Select which channels to copy to the new combined dataset'};

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''P''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'P'};

%--------------------------------------------------------------------------
% crop
%--------------------------------------------------------------------------
combineplanar          = cfg_exbranch;
combineplanar.tag      = 'combineplanar';
combineplanar.name     = 'Combine planar';
combineplanar.val      = {D, mode, prefix};
combineplanar.help     = {'Combine planar MEG channels'}';
combineplanar.prog     = @eeg_combineplanar;
combineplanar.vout     = @vout_eeg_combineplanar;
combineplanar.modality = {'EEG'};

%==========================================================================
function out = eeg_combineplanar(job)
% construct the S struct
S           = job;
S.D         = S.D{1};
out.D       = spm_eeg_combineplanar(S);
out.Dfname  = {fullfile(out.D)};

%==========================================================================
function dep = vout_eeg_combineplanar(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Planar-combined MEG data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Planar-combined MEG datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
