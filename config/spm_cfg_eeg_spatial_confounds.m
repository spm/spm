function sconfounds = spm_cfg_eeg_spatial_confounds
% configuration file for reading montage files
%__________________________________________________________________________
% Copyright (C) 2014-2016 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_spatial_confounds.m 6926 2016-11-09 22:13:19Z guillaume $

D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

timewin         = cfg_entry;
timewin.tag     = 'timewin';
timewin.name    = 'Time window';
timewin.strtype = 'r';
timewin.num     = [1 2];
timewin.val     = {[-Inf Inf]};
timewin.help    = {'Time window (ms)'};

condlabel = cfg_entry;
condlabel.tag = 'conditions';
condlabel.name = 'Condition label';
condlabel.strtype = 's';
condlabel.help = {''};

conditions = cfg_repeat;
conditions.tag = 'condrepeat';
conditions.name = 'Conditions';
conditions.help = {'Specify the labels of the conditions to include in the SVD.'};
conditions.num  = [0 Inf];
conditions.values  = {condlabel};
conditions.val = {};

ncomp = cfg_entry;
ncomp.tag = 'ncomp';
ncomp.name = 'Number of components';
ncomp.strtype = 'n';
ncomp.num = [1 1];
ncomp.help = {'Number of confound components to keep.'};

threshold = cfg_entry;
threshold.tag = 'threshold';
threshold.name = 'Threshold';
threshold.strtype = 'r';
threshold.val = {NaN};
threshold.num = [1 1];
threshold.help = {'Threshold for data amplitude after correction.',...
    'If defined this overrides number of components.'};

svd = cfg_branch;
svd.tag = 'svd';
svd.name = 'SVD';
svd.val = {timewin, conditions, threshold, ncomp};
svd.help = {'Define confounds from SVD of artefact samples.'};

conffile = cfg_files;
conffile.tag = 'conffile';
conffile.name = 'File name';
conffile.filter = 'mat';
conffile.num = [1 1];
conffile.help = {'Select the EEG mat file.'};

spmeeg = cfg_branch;
spmeeg.tag = 'spmeeg';
spmeeg.name = 'SPM M/EEG dataset';
spmeeg.val = {conffile};
spmeeg.help = {'Load confounds defined in a different dataset.'};

conffile = cfg_files;
conffile.tag = 'conffile';
conffile.name = 'BESA confounds file';
conffile.filter = '.*.bsa';
conffile.num = [1 1];
conffile.help = {'Select the EEG mat file.'};

besa = cfg_branch;
besa.tag = 'besa';
besa.name = 'BESA';
besa.val = {conffile};
besa.help = {'Load confounds defined in BESA from a .bsa file.'};

eyes = cfg_const;
eyes.tag = 'eyes';
eyes.name = 'Eyes';
eyes.val  = {1};
eyes.help = {'Defined confounds based on standard locations for the eyes',...
    'This is a very crude method, a data-driven method is preferable'};

clr = cfg_const;
clr.tag = 'clear';
clr.name = 'Clear';
clr.val  = {1};
clr.help = {'Clear previously defined spatial confounds.'};

mode = cfg_repeat;
mode.tag = 'mode';
mode.name = 'Mode';
mode.values = {svd, spmeeg, besa, eyes, clr};
mode.help = {'Select methods for defining spatial confounds.'};

sconfounds = cfg_exbranch;
sconfounds.tag = 'sconfounds';
sconfounds.name = 'Define spatial confounds';
sconfounds.val = {D, mode};
sconfounds.help = {'Define spatial confounds for topography-based correction of artefacts.'};
sconfounds.prog = @eeg_sconfounds;
sconfounds.vout = @vout_eeg_sconfounds;
sconfounds.modality = {'EEG'};

function out = eeg_sconfounds(job)

D = char(job.D{1});

for i = 1:numel(job.mode)
    mode = job.mode{i};
    S = [];
    S.D = D;
    S.mode = char(fieldnames(mode));
    switch S.mode
        case 'svd'
            if ~isnan(mode.svd.threshold)
                S.threshold = mode.svd.threshold;
            else
                S.ncomp = mode.svd.ncomp;
            end
            S.timewin = mode.svd.timewin;
        case {'besa', 'spmeeg'}
            S.conffile = char(mode.(S.mode).conffile);
        otherwise
            % do nothing
    end
    
    D = spm_eeg_spatial_confounds(S);
end


out.D = {fullfile(D)};

function dep = vout_eeg_sconfounds(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Dataset with spatial confounds';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});


