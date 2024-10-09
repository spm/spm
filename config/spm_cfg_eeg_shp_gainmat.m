function eeg_pca_gainmat = spm_cfg_eeg_shp_gainmat
% Configuration file for creating distorted versions of subject anatomy
% Based on original antomical and predetermined 100 eigen component template space.
%__________________________________________________________________________

% Gareth Barnes, Yael Balbastre
% Copyright (C) 2024 Imaging Neuroscience

eeg_pca_gainmat          = cfg_exbranch;
eeg_pca_gainmat.tag      = 'eeg_shp_gainmat';
eeg_pca_gainmat.name     = 'Gain matrices for surfaces';
eeg_pca_gainmat.val      = @eeg_shp_gainmat_cfg;
eeg_pca_gainmat.help     = {'To compute new lead field/ gain matrices for multiple disto'};
eeg_pca_gainmat.prog     = @specify_eeg_shp_gainmat;
eeg_pca_gainmat.vout     = @vout_specify_eeg_shp_gainmat;
eeg_pca_gainmat.modality = {'MEG'};


%==========================================================================
function varargout = eeg_shp_gainmat_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

D           = cfg_files;
D.tag       = 'D';
D.name      = 'M/EEG datasets';
D.filter    = 'mat';
D.val       = {''};
% D.dir       = {'C:\Users\gbarnes\Documents\jimmydata\output\data\gb070167\'};
D.num       = [1 1];
D.help      = {'Select the M/EEG mat file'};

val         = cfg_entry;
val.tag     = 'val';
val.name    = 'Inversion index';
val.strtype = 'n';
val.help    = {'Index of the cell in D.inv where the results will be stored.'};
val.val     = {1};


% surface_files           = cfg_files;
% surface_files.tag       = 'surface_files';
% surface_files.name      = 'Surface files ';
% surface_files.filter    = 'gii';
% surface_files.val       = {''};
% surface_files.dir       = {''};
% % surface_files.num       = [Inf Inf];
% surface_files.help      = {'Select the surface files (*.gii) for which to compute new lead fields'};

LFRedo                  = cfg_menu;
LFRedo.tag              = 'LFRedo';
LFRedo.name             = 'Force recompute of lead fields';
LFRedo.val              = {'No'};
LFRedo.help             = {'Will force recompute lead fields even if they exist'};
LFRedo.labels           = {'Yes', 'No'}';
LFRedo.values           = {'Yes', 'No'}';

LFheadmodel             = cfg_menu;
LFheadmodel.tag         = 'LFheadmodel';
LFheadmodel.name        = 'Use which head models';
LFheadmodel.val         = {'Single Shell'};
LFheadmodel.help        = {'Select the head model to use'};
LFheadmodel.labels      = {'Single Shell', 'Single Sphere'}';
LFheadmodel.values      = {'Single Shell', 'Single Sphere'}';


[cfg,varargout{1}] = deal({D, val,LFheadmodel,LFRedo});


%==========================================================================
function  out = specify_eeg_shp_gainmat(job)

out.D      = {};
val        = job.val;

% ------------------------------
% Load D and make a working copy
% ------------------------------
D          = spm_eeg_load(job.D{val});
folder_out = [D.path filesep 'outPCA' filesep];
Dc         = D.copy(folder_out); % work with copy of the original

% ------------------------
% Load all deformed meshes
% ------------------------
path_mri       = Dc.inv{1}.mesh.sMRI;
folder_mri     = spm_file(path_mri, 'path');
folder_mesh    = fullfile(folder_mri, 'PCA', 'Cerebros');
[brianmeshes]  = spm_select('FPList', folder_mesh, '.*\.gii$');
Nfiles         = size(brianmeshes,1);

fprintf('\n Found %d surfaces in %s \n',Nfiles,folder_mesh);
fprintf('\n Plus using original surface %s \n',Dc.inv{1}.mesh.tess_ctx);

% ------------------------------
% Prepare folder for head models
% ------------------------------
headmodel        = job.LFheadmodel;
% remove spaces for dir name
uheadmodel       = regexprep(headmodel,' ', '_'); 
% make a directory to contain headmodels and then surface leadfields
folder_headmodel = fullfile(folder_out, uheadmodel); 
mkdir(folder_headmodel);

% -------------------------------------------------------------------------
% Compute forward model from real brain
% -------------------------------------------------------------------------
% It actually does a bit more than that and sets up various fields 
% related to the head model and inversion model.

Dc.inv{val}.forward(1).voltype = headmodel;
Dc = spm_eeg_inv_forward(Dc);

% -------------------------------------------------------------------------
% Loop over true brain (0) and deformed brains ("brians", > 0)
% -------------------------------------------------------------------------
% Rather than recomputing everything using spm_eeg_inv_forward, we
% only update the cortical mesh and update the gain matrix

brainmesh = Dc.inv{val}.mesh.tess_ctx;
for cerebro = 0:Nfiles
    fprintf('\n Creating lead fields for brian %d of %d \n',cerebro, Nfiles);

    if cerebro == 0
        % Iteration 0 corresponds to the true brain !
        brianmesh = brainmesh;
    else
        % Distorted brain!
        % Keep file as it was on the first run
        % replace subject's cortex and delete any ref to lead fields
        brianmesh = deblank(brianmeshes(cerebro,:));
        Dc.inv{val}.mesh.tess_ctx = brianmesh;
    end

    meshname = spm_file(brianmesh, 'basename');

    % Load mesh (its vertices are in native space)
    gii_native = gifti(brianmesh);
    Vn = gii_native.vertices; 
    Nv = length(Vn);

    % Move from native MRI space to MEG space
    % (got this transform from spm_eeg_inv_forward)
    M    = Dc.inv{val}.datareg.fromMNI * Dc.inv{val}.mesh.Affine;
    M    = diag([1e-3 1e-3 1e-3 1]) * M;
    Vctf = (M*[gii_native.vertices ones(Nv,1)]')';

    % topology stays the same, so no need to update faces
    Dc.inv{val}.forward.mesh.vert = double(Vctf(:,1:3));

    % -----------------------------------------------------------
    % Run forward model from distored brain (only update gainmat)
    % -----------------------------------------------------------
    fprintf('\nComputing lead fields')

    path_gainmat = fullfile(folder_headmodel, ['GainMat_' meshname '.mat']);
    if ~exist(path_gainmat, 'file') || strcmp(job.LFRedo,'Yes')
        Dc.inv{val}.gainmat = '';
        Dc.save;
        % Generate/load lead field
        [~,Dc]= spm_eeg_lgainmat(Dc);			
        spm_copy(fullfile(Dc.path, Dc.inv{val}.gainmat),path_gainmat);
    end
    allgainmat{cerebro+1}=path_gainmat;
end


%--------------------------------------------------------------------------
% reset original mesh in the copied file
%--------------------------------------------------------------------------
Dc.inv{val}.forward.mesh.vert = Dc.inv{val}.forward.mesh.vert;
Dc.inv{val}.mesh.tess_ctx     = D.inv{val}.mesh.tess_ctx;
Dc.inv{val}.gainmat           = ''; % remove LF reference
Dc.save;

out.Dc{1}   = fullfile(Dc.path, Dc.fname); % output copy of original file
out.gainmat = allgainmat;


%==========================================================================
function dep = vout_specify_eeg_shp_gainmat(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) with a forward model';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
