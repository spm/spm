spm('defaults', 'eeg');

%% Basic definitions
% Change to the path to the multimodal MEG folder on your computer
root = 'D:\Data\Multimodal\CTF MEG';

% List of datasets (runs)
datasets = {'SPM_CTF_MEG_example_faces1_3D.ds', 'SPM_CTF_MEG_example_faces2_3D.ds'};

% List of condition labels,
condlabels  = {'faces', 'scrambled'};

% corresponding event types
eventtypes  = {'UPPT001', 'UPPT001'};

% and corresponding event values
eventvalues = {1, 2};

%% Read and pre-process the data
data = {};
condtrials = zeros(1, numel(condlabels));
for d = 1:numel(datasets)
    for c = 1:numel(condlabels)        
        
        cfg = [];
        
        % trial definition
        cfg.dataset = fullfile(root, datasets{d});
        cfg.trialdef.eventtype  = eventtypes{c};
        cfg.trialdef.eventvalue = eventvalues{c};
        cfg.trialdef.prestim    = 0.2; % from -200
        cfg.trialdef.poststim   = 0.6; % to 600 ms
        
        cfg = ft_definetrial(cfg);
        
        % adjust trigger latency (see the multimodal chapter)
        % here any other manipulation with the trl can be added
        hdr = ft_read_header(cfg.dataset);
        cfg.trl(:,1:2) = cfg.trl(:,1:2) + round(25*hdr.Fs/1000);
        
        cfg.channel = 'MEG';
        
        % baseline correction
        cfg.demean     = 'yes';
        cfg.baselinewindow = [-0.2 0];
        
        % actually reading the data
        data{d, c}  = ft_preprocessing(cfg);
        
        % remember how many trials for each condition
        condtrials(c) = condtrials(c) + numel(data{d, c}.trial);
    end
end
%% Put all the conditions and datasets in one struct 
% because the : operator take the data column-wise all the trials for the
% same condition will be together
data = ft_appenddata([], data{:});
%% Downsample
cfg = [];    
cfg.resamplefs = 200;
cfg.detrend    = 'no';
data = ft_resampledata(cfg, data);
%% Convert to SPM8 format
D = spm_eeg_ft2spm(data, ['ft_' spm_file(datasets{1}, 'basename')]);

%% Posp-processing of converted data

% Read sensors and fiducials from the first dataset
D = sensors(D, 'MEG', ft_convert_units(ft_read_sens(fullfile(root, datasets{1})), 'mm'));
D = fiducials(D, ft_convert_units(ft_read_headshape(fullfile(root, datasets{1})), 'mm'));

%% Set condition labels using the previously stored numbers of trials
condtrials = cumsum([1 condtrials]);

for c = 1:numel(condlabels)
    ind = condtrials(c):(condtrials(c+1)-1);
    if ~isempty(ind)
        D = conditions(D, ind, condlabels{c});
    end
end

%% Create 2D channel layout
S = [];
S.task = 'project3D';
S.modality = 'MEG';
S.updatehistory = 1;
S.D = D;
D = spm_eeg_prep(S);

%% Save the dataset
save(D);
%%
%% Create a head model
matlabbatch{1}.spm.meeg.source.headmodel.D = {fullfile(D.path, D.fname)};

matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 3;

matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'FIL_CTF_L';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'FIL_CTF_R';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;


matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';

spm_jobman('run', matlabbatch);
