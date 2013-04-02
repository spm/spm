function S = spm_cfg_eeg_contrast
% configuration file for computing contrast over epochs
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_contrast.m 5377 2013-04-02 17:07:57Z vladimir $

D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

c = cfg_entry;
c.tag = 'c';
c.name = 'Contrast coefficients';
c.strtype = 'r';
c.num = [1 inf];
c.help = {'Enter the contrast vector.'};

label = cfg_entry;
label.tag = 'label';
label.name = 'New condition label';
label.strtype = 's';
label.help = {'Enter the label for the condition derived by applying the contrast.'};

contrast = cfg_branch;
contrast.tag = 'contrast';
contrast.name = 'Contrast';
contrast.val = {c label};

contrasts = cfg_repeat;
contrasts.tag = 'contrasts';
contrasts.name = 'Contrasts';
contrasts.help = {'Each contrast defines a new condition in the output file.'};
contrasts.values  = {contrast};
contrasts.num     = [1 Inf];

weighted = cfg_menu;
weighted.tag = 'weighted';
weighted.name = 'Weight by replications';
weighted.labels = {'Yes', 'No'};
weighted.values = {1 , 0};
weighted.val = {1};
weighted.help = {'Weight the contrast by the numner of replications from which the average was computed.'};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''w''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'w'};

S = cfg_exbranch;
S.tag = 'contrast';
S.name = 'Contrast over epochs';
S.val = {D contrasts weighted prefix};
S.help = {'Computes contrasts over EEG/MEG epochs.'};
S.prog = @eeg_contrast;
S.vout = @vout_eeg_contrast;
S.modality = {'EEG'};

function out = eeg_contrast(job)
% construct the S struct
S.D = job.D{1};
S.c = cat(1, job.contrast(:).c);
S.label = {job.contrast.label};

S.weighted = job.weighted;
S.prefix   = job.prefix;

out.D = spm_eeg_contrast(S);

out.Dfname = {fullfile(out.D)};

function dep = vout_eeg_contrast(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Contrast of M/EEG epochs';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Contrast Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});

