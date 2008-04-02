function S = spm_cfg_eeg_artefact
% configuration file for M/EEG artefact detection
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_artefact.m 1295 2008-04-02 14:31:24Z volkmar $

rev = '$Rev: 1295 $';
D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

nothing = cfg_const;
nothing.tag = 'nothing';
nothing.name = 'No lists';
nothing.val{1} = double([]);

out_list = cfg_entry;
out_list.tag = 'out_list';
out_list.name = 'Artefactual trial list';
out_list.strtype = 'r';
out_list.num = [0 inf];
out_list.help = {'List artefactual trials (0 for none)'};

in_list = cfg_entry;
in_list.tag = 'in_list';
in_list.name = 'Clean trial list';
in_list.strtype = 'r';
in_list.num = [0 inf];
in_list.help = {'List clean trials (0 for none)'};

list = cfg_branch;
list.tag = 'list';
list.name = 'Trial lists';
list.val = {out_list in_list};

External_list         = cfg_choice;
External_list.tag     = 'External_list';
External_list.name    = 'Use own artefact list';
External_list.val = {nothing};
External_list.help    = {'Choose this option if you want to specify your own artefact list'}';
External_list.values = {list nothing};

Weightingfunction = cfg_entry;
Weightingfunction.tag = 'Weightingfunction';
Weightingfunction.name = 'Weightingfunction';
Weightingfunction.strtype = 'r';
Weightingfunction.num = [1 1];
Weightingfunction.help = {'Input offset of weighting function'};

Smoothing = cfg_entry;
Smoothing.tag = 'Smoothing';
Smoothing.name = 'Smoothing';
Smoothing.strtype = 'r';
Smoothing.num = [1 1];
Smoothing.help = {'FWHM for residual smoothing (ms)'};

weighted_arg = cfg_branch;
weighted_arg.tag = 'weighted_arg';
weighted_arg.name = 'Robust averaging';
weighted_arg.val = {Weightingfunction Smoothing};

nothing = cfg_const;
nothing.tag = 'nothing';
nothing.name = 'No robust averaging';
nothing.val{1} = double([]);

weighted         = cfg_choice;
weighted.tag     = 'weighted';
weighted.name    = 'Use robust averaging';
weighted.val = {nothing};
weighted.help    = {'Choose this option if you want to apply robust averaging'}';
weighted.values = {weighted_arg nothing};

nothing = cfg_const;
nothing.tag = 'nothing';
nothing.name = 'No thresholding';
nothing.val{1} = double([]);

threshold = cfg_entry;
threshold.tag = 'threshold';
threshold.name = 'Thresholds';
threshold.strtype = 'r';
threshold.num = [1 inf];
threshold.help = {'Channel-wise thresholds. Use single threshold to apply the same threshold to all channels'};

Check_Threshold         = cfg_choice;
Check_Threshold.tag     = 'Check_Threshold';
Check_Threshold.name    = 'Threshold channels';
Check_Threshold.val = {threshold};
Check_Threshold.help    = {'Choose this option if you want to threshold the data'}';
Check_Threshold.values = {threshold nothing};

artefact = cfg_branch;
artefact.tag = 'artefact';
artefact.name = 'Artefact';
artefact.val = {External_list weighted Check_Threshold};

S = cfg_exbranch;
S.tag = 'eeg_artefact';
S.name = 'M/EEG artefact detection';
S.val = {D artefact};
S.help = {'Detect artefacts in epoched M/EEG data.'};
S.prog = @eeg_artefact;
S.vout = @vout_eeg_artefact;
S.modality = {'EEG'};


function out = eeg_artefact(job)
% construct the S struct
S.D = job.D{1};
S.artefact = job.artefact;

if isfield(S.artefact.External_list, 'nothing')
    S.artefact.External_list = 0;
else
    S.artefact.External_list = 1;
    S.artefact.in_list = job.artefact.External_list.in_list;
    S.artefact.out_list = job.artefact.External_list.out_list;

end

if isfield(S.artefact.weighted, 'nothing')
    S.artefact.weighted = 0;
else
    S.artefact.weighted = 1;
    S.artefact.Weightingfunction = job.artefact.weighted.Weightingfunction;
    S.artefact.Smoothing = job.artefact.weighted.Smoothing;

end

if isfield(S.artefact.Check_Threshold, 'nothing')
    S.artefact.Check_Threshold = 0;
else
    S.artefact.Check_Threshold = 1;
    S.artefact.threshold = job.artefact.Check_Threshold.threshold;
end

out.D = spm_eeg_artefact(S);

function dep = vout_eeg_artefact(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Arteface detection';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});


