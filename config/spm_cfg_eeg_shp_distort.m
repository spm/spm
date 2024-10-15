function eeg_pca_distort = spm_cfg_eeg_shp_distort
% Configuration file for creating distorted versions of subject anatomy
% Based on original antomical and predetermined 100 eigen component template space.
%__________________________________________________________________________

% Gareth Barnes, Yael Balbastre
% Copyright (C) 2024 Imaging Neuroscience


eeg_pca_distort          = cfg_exbranch;
eeg_pca_distort.tag      = 'eeg_pca_distort';
eeg_pca_distort.name     = 'Distortion Trajectory';
eeg_pca_distort.val      = @eeg_shp_distort_cfg;
eeg_pca_distort.help     = {'To distort surface information in an M/EEG dataset based on typical normal variation'};
eeg_pca_distort.prog     = @specify_eeg_shp_distort;
eeg_pca_distort.vout     = @vout_specify_eeg_shp_distort;
eeg_pca_distort.modality = {'MEG'};


%==========================================================================
function varargout = eeg_shp_distort_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

D                       = cfg_files;
D.tag                   = 'D';
D.name                  = 'M/EEG datasets';
D.filter                = 'mat';
D.val                   = {''};
% D.dir                   = {'C:\Users\gbarnes\Documents\jimmydata\output\data\gb070167\'};
D.num                   =  [1 1];
D.help                  =  {'Select the M/EEG mat file'};

val                     = cfg_entry;
val.tag                 = 'val';
val.name                = 'Inversion index';
val.strtype             = 'n';
val.help                = {'Index of the cell in D.inv where the results will be stored.'};
val.val                 = {1};

PCAtemplate             = cfg_files;
PCAtemplate.tag         = 'PCAtemplate';
PCAtemplate.name        = 'PCA Template location';
PCAtemplate.filter      = '.*';
PCAtemplate.val         = {{fullfile(spm('Dir'),'tpm','shp','Template_0.nii')}};
% PCAtemplate.dir         = {fullfile(spm('Dir'),'tpm','shp')};
PCAtemplate.num         = [1 1];
PCAtemplate.help        = {'Select the generic PCA template file (Template_0.nii). Will download this automatically first time used'};

TemplateRedo            = cfg_menu;
TemplateRedo.tag        = 'TemplateRedo';
TemplateRedo.name       = 'Force recompute if exists';
TemplateRedo.val        = {'No'};
TemplateRedo.help       = {'Will force recompute template for this subject'}';
TemplateRedo.labels     = {'Yes', 'No'}';
TemplateRedo.values     = {'Yes', 'No'}';

DistortIndices          = cfg_entry;
DistortIndices.tag      = 'DistortIndices';
DistortIndices.name     = 'Range of indices';
DistortIndices.strtype  = 'i';
DistortIndices.num      = [1 2];
DistortIndices.val      = {[8 100]};
DistortIndices.help     = {'The range (from, to) of PCA indices to distort'};

Zcentre                 = cfg_menu;
Zcentre.tag             = 'Zcentre';
Zcentre.name            = 'Sample surfaces centered about the subject or about the population mean?';
Zcentre.val             = {'sub'};
Zcentre.help            = {
    'If sub(ject), the generated surfaces will have a latent code in Z[subject] +/- Zrange\n'
    'If can(onical), the generated surfaces will have a latent code in 0 +/- Zrange'
    }';
Zcentre.labels          = {'sub', 'can'}';
Zcentre.values          = {'sub', 'can'}';

Zrange                  = cfg_entry;
Zrange.tag              = 'Zrange';
Zrange.name             = 'Range of Z';
Zrange.strtype          = 'r';
Zrange.num              = [1 2];
Zrange.val              = {[0 3]};
Zrange.help             = {'The normal range (in Z) over which to perterb true brain (eg [0,+3])'};

Zlimit                 = cfg_entry;
Zlimit.tag             = 'Zlimit';
Zlimit.name            = 'Z max';
Zlimit.strtype         = 'r';
Zlimit.num             = [1 1];
Zlimit.val             = {3.5};
Zlimit.help            = {'The max absolute Z value permitted (will flip displacement sign otherwise)'};


Npoints                 = cfg_entry;
Npoints.tag             = 'Npoints';
Npoints.name            = 'Number points';
Npoints.strtype         = 'i';
Npoints.num             = [1 1];
Npoints.val             = {17};
Npoints.help            = {'The number of points on the distortion trajectory (linearly spaced between z limits)'};

RandSeed                = cfg_entry;
RandSeed.tag            = 'RandSeed';
RandSeed.name           = 'Random seed';
RandSeed.strtype        = 'i';
RandSeed.num            = [1 1];
RandSeed.val            = {1};
RandSeed.help           = {'The random seed that defines trajectory'};


[cfg,varargout{1}] = deal({D, val, PCAtemplate,TemplateRedo, DistortIndices, Npoints,Zrange,Zlimit,Zcentre,RandSeed});


%==========================================================================
function  out = specify_eeg_shp_distort(job)

out.D = {};

% Directory that contains shape model files
folder_shp = spm_file(job.PCAtemplate{1}, 'path');

% Number of points on distortion trajectory
K   = job.Npoints;
D   = spm_eeg_load(job.D{1});
val = job.val;

% Check everything is setup correctly
if ~isfield(D,'inv')
    error('no head model set up');
end

% Resets gainmat as these will be incorrect for the deformed anatomy
D.inv{val}.gainmat = '';

%D = spm_eeg_inv_forward(D);

% MRI in native space
path_mri   = D.inv{1}.mesh.sMRI;
path_mesh  = D.inv{1}.mesh.tess_ctx;
folder_mri = spm_file(path_mri, 'path');
name_mri   = spm_file(path_mri, 'basename');

% Folder where this subject's output PCA files are saved
% The space the MRI ends up in is not typical mni canonical space, but 
% the "shape template" space that was used to create the dictionary.
% Calling this "shp" space
folder_shp_out    = fullfile(folder_mri, 'PCA');
path_latent       = fullfile(folder_shp_out, ['pca_'    name_mri '.mat']);
path_velocity     = fullfile(folder_shp_out, ['iv_drc1' name_mri '_Template.nii']);
path_forward      = fullfile(folder_shp_out, [ 'y_drc1' name_mri '_Template.nii']);

% If redo, delete existing PCA directory
if strcmp(job.TemplateRedo,'Yes')
    if exist(folder_shp_out, 'dir')
        fprintf('\nRemoving existing PCA directory for fresh start\n')
        rmdir(folder_shp_out,'s');
    end
end

% If it hasn't been done yet, we must register the native mri with
% the shape template (Native -> Import -> Shape Model -> Latent code)
if ~exist(path_latent, 'file')
    fprintf('\nDid not find latent file, transforming\n')
    spm_shp_get_transforms(path_mri, [], folder_shp);
else
    fprintf('\nFound existing latent file for this subject not transforming\n')
end

idxPC = job.DistortIndices(1):job.DistortIndices(2);  % Components that will be deformed
rng(job.RandSeed, 'twister');                         % Same random deviations for everyone
sgnPC = sign(randn(1,length(idxPC)));                 % This sets distortion trajectory



% -------------------------------------------------------------------------
% Generate distorted brains
span = job.Zrange;
PC   = idxPC .* sgnPC;
fprintf(...
    '\nCreating a trajectory of %d brains from z=%3.2f to %3.2f (random signs) with seed %d\n', ...
    K, span(1), span(2), job.RandSeed ...
)
[z,z0,archivos] = spm_shp_sample_brains(         ...
    path_mesh,                                  ... Path to cortical mesh
    K,                                          ... Number of deformed surfaces 
    'can',    strcmpi(job.Zcentre, 'can'),      ... Deform about true brain
    'pc',     PC,                               ... Components to deform
    'span',   span,                             ... Deformation span (z scores)
    'zlimit', job.Zlimit,                           ... Max value of z allowed     
    'fout',   folder_shp_out,                   ... Output folder
    'fshp',   folder_shp,                       ... Model folder
    'suffix', sprintf('%03d', job.RandSeed),    ... Suffix (random seed)
    'z0',     path_latent,                      ... Subject's latent code
    'v0',     path_velocity,                    ... Subject's initial velocity
    'y0',     path_forward,                     ... Subject's forward transform
    'r2n',    path_latent                       ... Subject's import-to-native transform
);



% -------------------------------------------------------------------------
% Save and return
save(D);
out.D{1}   = fullfile(D.path, D.fname);
out.newz=z; %% new z scores of distortions (changed from PC to keep within z thresh)
out.z0=z0; %% latent code for original/ true brain.
out.brain  = D.inv{1}.mesh.tess_ctx;
out.brians = archivos;


%==========================================================================
function dep = vout_specify_eeg_shp_distort(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) with a forward model';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});