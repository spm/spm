function correct = spm_cfg_eeg_correct_sensor_data
% configuration file for reading montage files
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_correct_sensor_data.m 6029 2014-05-30 18:52:03Z vladimir $

D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};


mode = cfg_menu;
mode.tag = 'mode';
mode.name = 'Correction mode';
mode.labels = {'SSP', 'Berg'};
mode.values = {'ssp', 'berg'};
mode.val    = {'ssp'};
mode.help = {'Select correction method.',...
    'SSP removes more but also distorts the data more',...
    'Berg method requires forward model to be defined in the dataset'};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''T''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'T'};

correct = cfg_exbranch;
correct.tag = 'correct';
correct.name = 'Correct sensor data';
correct.val = {D, mode, prefix};
correct.help = {'Perform topography-based correction of artefacts'};
correct.prog = @eeg_correct;
correct.vout = @vout_eeg_correct;
correct.modality = {'EEG'};

function out = eeg_correct(job)
% construct the S struct

S   = job;
S.D = char(S.D);

D = spm_eeg_correct_sensor_data(S);
out.D = {fullfile(D)};

function dep = vout_eeg_correct(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Dataset with spatial confounds';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});


