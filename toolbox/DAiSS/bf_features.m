function out = bf_features
% Prepare data features for filter computation
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2015-2023 Wellcome Centre for Human Neuroimaging


out          = cfg_exbranch;
out.tag      = 'features';
out.name     = 'Covariance features';
out.val      = @bf_features_cfg;
out.help     = {'Define features for covariance computation'};
out.prog     = @bf_features_run;
out.vout     = @bf_features_vout;
out.modality = {'EEG'};


%==========================================================================
function varargout = bf_features_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

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

cross_terms = cfg_menu;
cross_terms.tag = 'cross_terms';
cross_terms.name = 'Zero cross terms';
cross_terms.help = {'When fusing, set cross-terms between modalities to zero'};
cross_terms.labels = {'MEG to EEG only' 'MEG, MEGPLANAR, EEG', 'No'};
cross_terms.values = {'megeeg', 'all', 'no'};
cross_terms.val = {'megeeg'};

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

visualise = cfg_menu;
visualise.tag = 'visualise';
visualise.name = 'Visualise eigen-spectrum';
visualise.help = {'Visualise covariance log eigen-spectrum to check for effective data dimensionality.'};
visualise.labels = {'yes', 'no'};
visualise.values = {1, 0};
visualise.val = {1};

[cfg,varargout{1}] = deal({BF, whatconditions, woi, modality, fuse, cross_terms, plugin, reg, bootstrap, visualise});


%==========================================================================
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
end

if isempty(classchanind)
    %%%%%%%%%%%%
    S(1).samples = {};
    
    for i = 1:size(job.woi, 1)
        S.samples{i} = D.indsample(1e-3*job.woi(i, 1)):D.indsample(1e-3*job.woi(i, 2));
        if isnan(S.samples{i})
            error('Window specified not in dataset');
        end
    end
    %%%%%%%%%%%%
    % MWW 19/11/2014
else
    try
        classchanind=find(strcmp(D.chanlabels,'Class')); % MWW 19/11/2014
    catch
        error('There must be a Class channel in D if job.woi is not specfied');
    end
end
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
    if ~isa(modalities{m}, 'cell')
        cmod = modalities(m);
    else
        cmod = modalities{m};
    end
    
    eegind = strmatch('EEG', cmod, 'exact');
    megind = [strmatch('MEG', cmod, 'exact') strmatch('MEGMAG', cmod, 'exact')];
    planarind = strmatch('MEGPLANAR', cmod, 'exact');
    
    chanind = cell(1, numel(cmod));
    chanind_cov = cell(1, numel(cmod));
    for n = 1:numel(cmod)
        chanind{n} = indchantype(BF.data.D, cmod{n}, 'GOOD');
        chanind_cov{n} = 1:length(chanind{n});
        if n>1
            chanind_cov{n} = chanind_cov{n} + length([chanind{1:(n-1)}]);
        end
    end
    
    S.channels=[chanind{:}];
    
    if isempty(S.channels)
        if ~isa(modalities{m}, 'cell')
            error('No good %s channels were found.\n', cmod{:});
        end
    end
    
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
        
        if numel(cmod)>1
            if isequal(job.cross_terms, 'all') && ~isempty(megind) && ~isempty(planarind)
                BF.features.(modality_name).C(chanind_cov{megind}, chanind_cov{planarind}) = 0;
                BF.features.(modality_name).C(chanind_cov{planarind}, chanind_cov{megind}) = 0;
            end
            if isequal(job.cross_terms, 'megeeg') && ~isempty(megind) && ~isempty(eegind)
                BF.features.(modality_name).C(chanind_cov{megind}, chanind_cov{eegind}) = 0;
                BF.features.(modality_name).C(chanind_cov{eegind}, chanind_cov{megind}) = 0;
            end
            if isequal(job.cross_terms, 'megeeg') && ~isempty(planarind) && ~isempty(eegind)
                BF.features.(modality_name).C(chanind_cov{planarind}, chanind_cov{eegind}) = 0;
                BF.features.(modality_name).C(chanind_cov{eegind}, chanind_cov{planarind}) = 0;
            end
        end
        
        
        if job.visualise
            F = spm_figure('GetWin', ['Log-eigenspectrum for ' modality_name]);clf;
            
            [~, SS, ~] = svd(BF.features.(modality_name).C);
            
            semilogy(diag(SS), '-o');
            
            xlabel('Singular value index');
            ylabel('Singular value');
        end
        
        S1.modality = modality_name;
        S1.chanind  = S.channels;
        
        BF.features.(modality_name) = feval(['bf_regularise_' reg_name], BF, S1);
        
        if job.visualise
            switch reg_name
                case {'clifftrunc','minkatrunc','mantrunc'}
                    switch reg_name
                        case {'clifftrunc','minkatrunc'}
                            plt = 1;
                            dof = size(BF.features.(modality_name).C,1);
                            str = sprintf('Estimated Rank: %03d',dof);
                        case {'mantrunc'}
                            plt = 1;
                            dof = size(BF.features.(modality_name).C,1);
                            str = sprintf('Specified Rank: %03d',dof);
                    end
                    if plt
                        if ~isempty(which('xline'))
                            xline(dof,'--',str);
                        else
                            try
                                hx = graph2d.constantline(dof, 'LineStyle',':', 'Color',[0 0 0]);
                                changedependvar(hx,'x');
                            end
                        end
                    end
                case 'manual'
                    [~, SS, ~] = svd(BF.features.(modality_name).C);
                    hold on
                    semilogy(diag(SS), '-o');
                    legend('Orignal','Regularised')
            end
        end
        
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
        for ii=1:NK
            
            S.samples = (D(classchanind,:,:)==ii);
            
            BF.data.samples.(modality_name).class{ii}=S.samples;
            
            BF.features.(modality_name).class{ii} = feval(['bf_features_' plugin_name], BF, S);
            
            
            if numel(cmod)>1
                if isequal(job.cross_terms, 'all') && ~isempty(megind) && ~isempty(planarind)
                    BF.features.(modality_name).class{ii}.C(chanind_cov{megind}, chanind_cov{planarind}) = 0;
                    BF.features.(modality_name).class{ii}.C(chanind_cov{planarind}, chanind_cov{megind}) = 0;
                end
                if isequal(job.cross_terms, 'megeeg') && ~isempty(megind) && ~isempty(eegind)
                    BF.features.(modality_name).class{ii}.C(chanind_cov{megind}, chanind_cov{eegind}) = 0;
                    BF.features.(modality_name).class{ii}.C(chanind_cov{eegind}, chanind_cov{megind}) = 0;
                end
                if isequal(job.cross_terms, 'megeeg') && ~isempty(planarind) && ~isempty(eegind)
                    BF.features.(modality_name).class{ii}.C(chanind_cov{planarind}, chanind_cov{eegind}) = 0;
                    BF.features.(modality_name).class{ii}.C(chanind_cov{eegind}, chanind_cov{planarind}) = 0;
                end
            end
            
            S1.chanind  = S.channels;
            S1.modality=modality_name;
            S1.class=ii;
            BF.features.(modality_name).class{ii} = feval(['bf_regularise_' reg_name], BF, S1);
            
        end
        
    end
    %%%%%%%%%%%%
    
    BF.features.(modality_name).chanind = S.channels;
end

BF.features.trials = S.trials;

bf_save_path(BF,fullfile(outdir, 'BF.mat'));

out.BF{1} = fullfile(outdir, 'BF.mat');


%==========================================================================
function dep = bf_features_vout(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
