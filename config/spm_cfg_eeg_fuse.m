function fuse = spm_cfg_eeg_fuse
% configuration file for fusing M/EEG files
%_______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_fuse.m 5377 2013-04-02 17:07:57Z vladimir $

rev = '$Rev: 5377 $';
D = cfg_files;
D.tag = 'D';
D.name = 'File Names';
D.filter = 'mat';
D.num = [2 Inf];
D.help = {'Select the M/EEG mat files.'};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the fused dataset. Default prefix is ''u''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'u'};

fuse = cfg_exbranch;
fuse.tag = 'fuse';
fuse.name = 'Fusion';
fuse.val = {D, prefix};
fuse.help = {'Fuse EEG/MEG data.'};
fuse.prog = @eeg_fuse;
fuse.vout = @vout_eeg_fuse;
fuse.modality = {'EEG'};

function out = eeg_fuse(job)
% construct the S struct
S = job;
S.D = char(S.D);

out.D = spm_eeg_fuse(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

function dep = vout_eeg_fuse(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Fused Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Fused Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});


