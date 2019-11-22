function out = bf_features
% Prepares data features for filter computation
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_features.m 7703 2019-11-22 12:06:29Z guillaume $

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
woi.name = 'Time windows of interest';
woi.strtype = 'r';
woi.num = [Inf 2];
woi.val = {[-Inf Inf]};
woi.help = {'Time windows to average over (ms)'};

modality = cfg_menu;
modality.tag = 'modality';
modality.name = 'Select modalities';
modality.help = {'Select modalities for the inversion (only relevant for multimodal datasets).'};
modality.labels = {'All', 'EEG', 'MEG', 'MEGPLANAR', 'EEG+MEG', 'MEG+MEGPLANAR', 'EEG+MEGPLANAR'};
modality.values = {
    {'EEG', 'MEG', 'MEGPLANAR'}
    {'EEG'}
    {'MEG'}
    {'MEGPLANAR'}
    {'EEG', 'MEG'}
    {'MEG', 'MEGPLANAR'}
    {'EEG', 'MEGPLANAR'}
    }';
modality.val = {{'MEG'}};


fuse = cfg_menu;
fuse.tag = 'fuse';
fuse.name = 'Fuse modalities';
fuse.help = {'Fuse sensors for different modalities together (requires prior rescaling).'};
fuse.labels = {'Don''t fuse' 'Fuse MEG only', 'Fuse all'};
fuse.values = {'no', 'meg', 'all'};
fuse.val = {'no'};

%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
plugin      = cfg_choice;
plugin.tag  = 'plugin';
plugin.name = 'Covariance computation method';

feature_funs = spm_select('List', fileparts(mfilename('fullpath')), '^bf_features_.*\.m$');
feature_funs = cellstr(feature_funs );
for i = 1:numel(feature_funs)
    plugin.values{i} = feval(spm_file(feature_funs{i},'basename'));
end


%--------------------------------------------------------------------------
% regularisation/reduction
%--------------------------------------------------------------------------
reg         = cfg_choice;
reg.tag  = 'regularisation';
reg.name = 'Regularisation method';

reg_funs = spm_select('List', fileparts(mfilename('fullpath')), '^bf_regularise_.*\.m$');
reg_funs = cellstr(reg_funs );
for i = 1:numel(reg_funs)
    reg.values{i} = feval(spm_file(reg_funs{i},'basename'));
end

bootstrap = cfg_menu;
bootstrap.tag = 'bootstrap';
bootstrap.name = 'Bootstrap';
bootstrap.labels = {'yes', 'no'};
bootstrap.values = {true, false};
bootstrap.val = {false};

out = cfg_exbranch;
out.tag = 'features';
out.name = 'Covariance features';
out.val = {BF, whatconditions, woi, modality, fuse, plugin, reg, bootstrap};
out.help = {'Define features for covariance computation'};
out.prog = @bf_features_run;
out.vout = @bf_features_vout;
out.modality = {'EEG'};
end

function  out = bf_features_run(job)

outdir = spm_file(job.BF{1}, 'fpath');

% cd(outdir);

BF = bf_load(fullfile(outdir, 'BF.mat'));
D  = BF.data.D;


plugin_name = cell2mat(fieldnames(job.plugin));
S         = job.plugin.(plugin_name);

%%%%%%%%%%%%
% MWW 19/11/2014   
classchanind=[];
try
    classchanind=find(strcmp(D.chanlabels,'Class')); % MWW 19/11/2014
catch
end;

if isempty(classchanind)
%%%%%%%%%%%%
    S(1).samples = {};

    for i = 1:size(job.woi, 1)
        S.samples{i} = D.indsample(1e-3*job.woi(i, 1)):D.indsample(1e-3*job.woi(i, 2));
        if isnan(S.samples{i})
            error('Window specified not in dataset');
        end;
    end
%%%%%%%%%%%%
% MWW 19/11/2014                
else
    try
        classchanind=find(strcmp(D.chanlabels,'Class')); % MWW 19/11/2014
    catch
        error('There must be a Class channel in D if job.woi is not specfied'); 
    end; 
end;
%%%%%%%%%%%%

if isfield(job.whatconditions, 'all')
    S(1).trials = D.indtrial(D.condlist, 'GOOD');
else    
    S(1).trials = D.indtrial(job.whatconditions.condlabel, 'GOOD');
    if isempty(S.trials)
        error('No trials matched the selection, check the specified condition labels');
    end
end

if job.bootstrap
    S.trials = S.trials(ceil(rand(1, length(S.trials)).*length(S.trials)));
end

reg_name   = cell2mat(fieldnames(job.regularisation));
S1         = job.regularisation.(reg_name);

switch job.fuse
    case 'no'
        modalities = job.modality;
    case 'meg'
        modalities{1} = intersect(job.modality, {'MEG', 'MEGMAG', 'MEGPLANAR'});%% added MEGMAG
        modalities    = [modalities intersect(job.modality, {'EEG'})];
    case 'all'
        modalities{1} = job.modality;
end

for m = 1:numel(modalities)
    chanind = indchantype(BF.data.D, modalities{m}, 'GOOD');
    if isempty(chanind)
        error(['No good ' modalities{m} ' channels were found.']);
    end
    S.channels=chanind;
    
    if isequal(char(modalities{m}), 'EEG') 
        modality_name  = 'EEG';
    elseif isequal(char(modalities{m}), 'MEGPLANAR')
        modality_name  = 'MEGPLANAR';
    elseif isequal(char(modalities{m}), 'MEGMAG') 
        modality_name  = 'MEGMAG';    
    else
        modality_name  = 'MEG';
    end
    
    %%%%%%%%%%%%
    % MWW 19/11/2014                
    if isempty(classchanind)       
    %%%%%%%%%%%%
    
        BF.features.(modality_name) = feval(['bf_features_' plugin_name], BF, S);

        S1.modality = modality_name;
        S1.chanind  = chanind;

        BF.features.(modality_name) = feval(['bf_regularise_' reg_name], BF, S1);

    %%%%%%%%%%%%
    % MWW 19/11/2014
    % added to allow S.samples to be specified via a "Class" channel, which
    % specifies which time points correspond to each class. This is so that
    % class-specific features can be calculated separately using just the 
    % timepoints for each class. This can then be used later for doing 
    % source reconstruction specific to each class.   
    % For example, the class specific features (covariances) can be used
    % by bf_features_cov_bysamples and bf_inverse_lcmv_multicov
    else

        disp('Ignoring job.woi. Using Class channel in D object to determine the time samples to use');
        NK=max(squash(D(classchanind,:,:)));
        for ii=1:NK,

            S.samples = (D(classchanind,:,:)==ii);

            BF.data.samples.(modality_name).class{ii}=S.samples;

            BF.features.(modality_name).class{ii} = feval(['bf_features_' plugin_name], BF, S);

            S1.modality=modality_name;
            S1.class=ii;
            BF.features.(modality_name).class{ii} = feval(['bf_regularise_' reg_name], BF, S1);

        end;   
        
    end; 
    %%%%%%%%%%%%
    
    BF.features.(modality_name).chanind = chanind;
end

BF.features.trials = S.trials;

bf_save_path(BF,fullfile(outdir, 'BF.mat'));

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
