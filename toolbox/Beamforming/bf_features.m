function out = bf_features
% Prepares data features for filter computation
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_features.m 4847 2012-08-16 17:29:23Z vladimir $

% dir Directory
% ---------------------------------------------------------------------
BF = cfg_files;
BF.tag = 'BF';
BF.name = 'BF.mat file';
BF.filter = '^BF.mat$';
BF.num = [1 1];
BF.help = {'Select BF.mat file.'};

all = cfg_const;
all.tag = 'all';
all.name = 'All';
all.val  = {1};

condlabel = cfg_entry;
condlabel.tag = 'condlabel';
condlabel.name = 'Condition label';
condlabel.strtype = 's';
condlabel.val = {''};

conditions = cfg_repeat;
conditions.tag = 'conditions';
conditions.name = 'Conditions';
conditions.help = {'Specify the labels of the conditions to be included in the inversion'};
conditions.num  = [1 Inf];
conditions.values  = {condlabel};
conditions.val = {condlabel};

whatconditions = cfg_choice;
whatconditions.tag = 'whatconditions';
whatconditions.name = 'What conditions to include?';
whatconditions.values = {all, conditions};
whatconditions.val = {all};

woi = cfg_entry;
woi.tag = 'woi';
woi.name = 'Time window of interest';
woi.strtype = 'r';
woi.num = [Inf 2];
woi.val = {[-Inf Inf]};
woi.help = {'Time window to average over (ms)'};

foi = cfg_entry;
foi.tag = 'foi';
foi.name = 'Frequency window of interest';
foi.strtype = 'r';
foi.num = [1 2];
foi.val = {[0 0]};
foi.help = {'Frequency window (Hz)'};

%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
plugin      = cfg_choice;
plugin.tag  = 'plugin';
plugin.name = 'Covariance computation method';

feature_funs = spm_select('List', fullfile(spm('dir'), 'toolbox', 'Beamforming'), '^bf_features_.*\.m$');
feature_funs = cellstr(feature_funs );
for i = 1:numel(feature_funs)
    plugin.values{i} = feval(spm_file(feature_funs{i},'basename'));
end

out = cfg_exbranch;
out.tag = 'bf_features';
out.name = 'Covariance features';
out.val = {BF, whatconditions, woi, foi, plugin};
out.help = {'Define features for covariance computation'};
out.prog = @bf_features_run;
out.vout = @bf_features_vout;
out.modality = {'EEG'};
end

function  out = bf_features_run(job)

outdir = spm_file(job.BF{1}, 'fpath');

cd(outdir);

BF = bf_load('BF.mat', {'data'});
D  = BF.data.D;

plugin_name = cell2mat(fieldnames(job.plugin));

S         = job.plugin.(plugin_name);
S.samples = {};
for i = 1:size(job.woi, 1)
    S.samples{i} = D.indsample(job.woi(i, 1)):D.indsample(job.woi(i, 2));
end

if isfield(job.whatconditions, 'all')
    S.trials = 1:D.ntrials;
else    
    S.trials = D.pickconditions(job.whatconditions.condlabel);
    if isempty(S.trials)
        error('No trials matched the selectsion, check the specified condition labels');
    end
end

S.foi = job.foi;

modalities = {'MEG', 'EEG'};

for m = 1:numel(modalities)
    
    if isfield(BF.data, modalities{m})
        S.channels = setdiff(meegchannels(D, modalities{m}), D.badchannels);
        
        BF.features.C.(modalities{m}) = feval(['bf_features_' plugin_name], BF, S);
    end
end

bf_save(BF);

out.BF{1} = fullfile(outdir, 'BF.mat');
end

function dep = bf_features_vout(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
end
