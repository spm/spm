function [z,folder_out] = spm_shp_get_transforms(path_mri,folder_out,folder_shp)
% Transform MRI to model space and compute its latent code
%
% FORMAT spm_shp_get_transforms(path_mri,[folder_out],[folder_shp])
%
% path_mri   - Path to input MRI
% folder_out - Path to output folder (default: {input_folder}/PCA)
% folder_shp - Path to template folder (default: {spm('Dir')}/tpm/shp)
% z          - Latent vectors describing the space in which the MRI lives
% 
% The following files are written in folder_out:
% * pca_{mri_name}.mat - Shape representation of the MRI
%                           'z'   - Latent vector
%                           'r2n' - Import to native affine transform
%                           'n2r' - Native to import affine transform
% *    {mri_name}.nii  - Copy of input MRI
% *  v_{mri_name}.nii  - Initial velocity (in import space)
% * iv_{mri_name}.nii  - Initial velocity (in group space)
% *  y_{mri_name}.nii  - Import to group nonlinear transform
% * iy_{mri_name}.nii  - Group to import nonlinear transform
%__________________________________________________________________________

% Gareth Barnes, Yael Balbastre
% Copyright (C) 2024 Wellcome Centre for Human Neuroimaging

folder_tpm = fullfile(spm('Dir'), 'tpm');           % SPM Template
folder_inp = spm_file(path_mri, 'path');            % Input folder
folder_shp = spm_shp_install(folder_shp);           % Model folder
if nargin < 2 || isempty(folder_out)
    folder_out = fullfile(folder_inp, 'PCA');       % Output folder
end

%% take native MRI and make a copy

% remove volume index
path_mri = spm_file(path_mri, 'ext', spm_file(path_mri, 'ext'));

mkdir(folder_out);
spm_copy(path_mri, folder_out, 'nifti', true);
path_mri = spm_file(path_mri, 'path', folder_out);
mri_name = spm_file(path_mri, 'basename');

%% Append file separator to make batch nicer to read

folder_tpm = [folder_tpm filesep];
folder_shp = [folder_shp filesep];
folder_out = [folder_out filesep];

%% (1) Segment
% Segment the registered MRI
% made sure to output the native (c*) and normalized (rc*) segmentations

spm_jobman('initcfg')

preproc.channel.vols        = {[path_mri ',1']};                % Data volumes
preproc.channel.biasreg     = 0.001;	                        % Bias regularisation
preproc.channel.biasfwhm    = 60;	                            % Bias FWHM: 60mm cutoff
preproc.channel.write       = [1 1];	                        % Save nothing
preproc.tissue(1).tpm       = {[folder_tpm 'TPM.nii,1']};	    % Tissue probability map
preproc.tissue(1).ngaus     = 1;			                    % Num of Gaussians
preproc.tissue(1).native    = [1 1];	                        % Native tissue: Native + Dartel
preproc.tissue(1).warped    = [0 0];	                        % Save none
preproc.tissue(2).tpm       = {[folder_tpm 'TPM.nii,2']};
preproc.tissue(2).ngaus     = 1;
preproc.tissue(2).native    = [1 1];
preproc.tissue(2).warped    = [0 0];
preproc.tissue(3).tpm       = {[folder_tpm 'TPM.nii,3']};
preproc.tissue(3).ngaus     = 2;
preproc.tissue(3).native    = [1 1];
preproc.tissue(3).warped    = [0 0];
preproc.tissue(4).tpm       = {[folder_tpm 'TPM.nii,4']};
preproc.tissue(4).ngaus     = 3;
preproc.tissue(4).native    = [1 1];
preproc.tissue(4).warped    = [0 0];
preproc.tissue(5).tpm       = {[folder_tpm 'TPM.nii,5']};
preproc.tissue(5).ngaus     = 4;
preproc.tissue(5).native    = [1 1];
preproc.tissue(5).warped    = [0 0];
preproc.tissue(6).tpm       = {[folder_tpm 'TPM.nii,6']};
preproc.tissue(6).ngaus     = 2;
preproc.tissue(6).native    = [0 0];
preproc.tissue(6).warped    = [0 0];
preproc.warp.mrf            = 1;
preproc.warp.cleanup        = 1;
preproc.warp.reg            = [0 0.001 0.5 0.05 0.2];	    % Warping regularisation
preproc.warp.affreg         = 'mni';		                % Affine regularisation
preproc.warp.fwhm           = 0;			                % Smoothness
preproc.warp.samp           = 3;			                % Sampling distance
preproc.warp.write          = [0 0];		                % Deformation fields
preproc.warp.vox            = NaN;
preproc.warp.bb             = [NaN NaN NaN; NaN NaN NaN];

matlabbatch{1}.spm.spatial.preproc = preproc;

%% (2) Downsample
% Our shape model was built at a slightly lower resolution (2mm) than
% typical "dartel import" images (1.5 mm), so we need to further
% resample the latter.

coreg_write.ref    = {[folder_shp 'Template_4.nii,1']};
coreg_write.source = {
    [folder_out 'rc1' mri_name '.nii,1']
    [folder_out 'rc2' mri_name '.nii,1']
};
coreg_write.roptions.interp = 4;
coreg_write.roptions.wrap   = [0 0 0];
coreg_write.roptions.mask   = 0;
coreg_write.roptions.prefix = 'd';

matlabbatch{2}.spm.spatial.coreg.write = coreg_write;

%% (3) Shoot
% Use the Shoot toolbox in SPM to register drc* to Template*
% This will output the v* (velocity) y* (exponentiated transform) j* (jacobian) files

shoot.warp1.images = {
    {[folder_out 'drc1' mri_name '.nii,1']}
    {[folder_out 'drc2' mri_name '.nii,1']}
}';
shoot.warp1.templates = {
    [folder_shp 'Template_0.nii']
    [folder_shp 'Template_1.nii']
    [folder_shp 'Template_2.nii']
    [folder_shp 'Template_3.nii']
    [folder_shp 'Template_4.nii']
};

matlabbatch{3}.spm.tools.shoot = shoot;

spm('defaults', 'PET');
spm_jobman('run', matlabbatch);

%% Affine map from native to import space

fprintf('Compute Native-to-Import tranform\n');

nii  = nifti([folder_out 'rc1' mri_name '.nii']);
mat  = nii.mat;   % Import space
mat0 = nii.mat0;  % Native space

n2r  = mat/mat0;
r2n  = mat0/mat;

%% Delete unused files

fprintf('Delete unused files\n');

delete([folder_out 'BiasField_' mri_name '.nii']);
delete([folder_out              mri_name '_seg8.mat']);
delete([folder_out          'm' mri_name '.nii']);
for i1=1:5
    delete([folder_out 'c'   num2str(i1) mri_name '.nii']);
    delete([folder_out 'rc'  num2str(i1) mri_name '.nii']);
end
delete([folder_out    'drc1' mri_name '.nii']);
delete([folder_out    'drc2' mri_name '.nii']);
delete([folder_out 'j_drc1'  mri_name '_Template.nii']);

%% We need iv* (inverse velocity) and iy* (inverse transform)

fprintf('Shoot inverse\n');

path_v_dr2tpl = fullfile(folder_out, [ 'v_drc1' mri_name '_Template.nii']);
path_v_tpl2dr = fullfile(folder_out, ['iv_drc1' mri_name '_Template.nii']);
path_y_dr2tpl = fullfile(folder_out, [ 'y_drc1' mri_name '_Template.nii']);
path_y_tpl2dr = fullfile(folder_out, ['iy_drc1' mri_name '_Template.nii']);

% Load velocity (and drop time dimension)
nii_v_dr2tpl         = nifti(path_v_dr2tpl);
fullsize             = nii_v_dr2tpl.dat.dim;
nii_v_dr2tpl.dat.dim = nii_v_dr2tpl.dat.dim([1 2 3 5]);

% Geodesic Shooting
dft         = spm_shoot_defaults;
v           = single(nii_v_dr2tpl.dat());
[y,~,iv,iy] = spm_shoot3d(v, [2 2 2 dft.rparam], Inf);
% iv          = -iv;

% Save inverse velocity (and add back time dimension)
nii_v_tpl2dr                = nii_v_dr2tpl;
nii_v_tpl2dr.dat.dim        = fullsize;
nii_v_tpl2dr.dat.fname      = path_v_tpl2dr;
create(nii_v_tpl2dr);
nii_v_tpl2dr.dat(:,:,:,:,:) = reshape(iv, fullsize);

% Model orientation
mat = nii_v_tpl2dr.mat;

% Save forward transform
y                           = spm_shp_warps('compose', mat, y);
nii_y_dr2tpl                = nifti(path_v_dr2tpl);
nii_y_dr2tpl.dat.fname      = path_y_dr2tpl;
create(nii_y_dr2tpl);
nii_y_dr2tpl.dat(:,:,:,:,:) = reshape(y, fullsize);

% Save inverse transform
iy                          = spm_shp_warps('compose', mat, iy);
nii_y_tpl2dr                = nifti(path_v_tpl2dr);
nii_y_tpl2dr.dat.fname      = path_y_tpl2dr;
create(nii_y_tpl2dr);
nii_y_tpl2dr.dat(:,:,:,:,:) = reshape(iy, fullsize);

%% Latent code

fprintf('Compute latent code\n');

% Load inverse velocity
path_iv_dr2tpl        = fullfile(folder_out, ['iv_drc1' mri_name '_Template.nii']);
nii_iv_tpl2dr         = nifti(path_iv_dr2tpl);
nii_iv_tpl2dr.dat.dim = nii_iv_tpl2dr.dat.dim([1 2 3 5]);
iv_tpl2dr             = single(nii_iv_tpl2dr.dat());

% Compute latent code
fsubspace	= fullfile(folder_shp, 'subspace_scaled.nii');
fmodel		= fullfile(folder_shp, 'model_variables.mat');
z           = spm_shp_project_velocity(iv_tpl2dr,fmodel,fsubspace);

%% Save to disk

fprintf('Save latent code and affine transforms\n');

save(fullfile(folder_out, ['pca_' mri_name '.mat']), 'z', 'n2r', 'r2n');

%% Delete temporary files
% delete([folder_out 'iv_drc1' mri_name '_Template.nii']);
% delete([folder_out 'v_drc1'  mri_name '_Template.nii']);




