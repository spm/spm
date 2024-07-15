function [z,output_folder]=spm_pca_get_transforms(nativeMRIname,templates_folder,output_folder)
%% z are latent vectors describing the space in which the MRI lives : vector of 100 elements
%% also get a copy of nativeMRI in template space in output_folder (optional)

if nargin<3
    output_folder=[];
end;


[pwd1,subname,c1]=fileparts(nativeMRIname);

if isempty(output_folder),
    output_folder=[pwd1 filesep 'PCA' filesep];
end;


MRIname=[subname c1];



% Manually create 'original' folder with the MRI and Template folder

templates_folder=[templates_folder filesep];
tpl		= [templates_folder 'Template_4.nii'];				% Template brain

mkdir(output_folder);
orig_folder		= [output_folder 'original'];
fprintf('\n Making copy of original MRI in directory %s\n',orig_folder)
mkdir(orig_folder);


TPM_folder		= fullfile(spm('Dir'),['tpm' filesep]);
CAN_folder=fullfile(spm('Dir'),['canonical' filesep]);


%% take native MRIand make a copy

[a1,subname,c1]=fileparts(nativeMRIname);
copyfile([a1 filesep subname  '.*'],orig_folder);%% of original native MRI

subj_file=[subname c1];
copyfile([a1 filesep subname  '.*'],output_folder); %% put a copy in output folder also (will be modified).

% Rename to match Yael coding


%% Coregister
% Register the native MRI with the new (database specific) SPM template (TPM.nii)
%% this updates the header of the MRI in the output directory so that it now in the new TPM template space
%% (close to canonical spm space, but not quite)

spm_jobman('initcfg')

matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[TPM_folder 'TPM.nii,1']};	% Reference image
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[output_folder subj_file ',1']};	% Source image
    %'C:\Users\gbarnes\Documents\jimmydata\data\mri\gb070167\pd_nw_mtflash3d_v3e_HeadCast_12ch_0002\PCA\mmsMQ0484_orig.img,1,1'
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi'; % Objective function: Normalised mutual information
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2]; % Separation
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]; % Tolerances
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];	% Histogram smoothing

spm('defaults', 'PET');
spm_jobman('run', matlabbatch);
clear matlabbatch


%% Segment
% Segment the registered MRI
% made sure to output the native (c*) and normalized (rc*) segmentations

spm_jobman('initcfg')

matlabbatch{1}.spm.spatial.preproc.channel.vols = {[output_folder subj_file ',1']};	% Data volumes
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;	% Bias regularisation
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;	% Bias FWHM: 60mm cutoff
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];	% Save nothing
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[TPM_folder 'TPM.nii,1']};	% Tissue probability map
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;			% Num of Gaussins
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];	% Native tissue: Native + Dartel
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];	% Save none
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[TPM_folder 'TPM.nii,2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[TPM_folder 'TPM.nii,3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[TPM_folder 'TPM.nii,4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[TPM_folder 'TPM.nii,5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[TPM_folder 'TPM.nii,6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];	% Warping regularisation
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';		% Affine regularisation
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;			% Smoothness
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;			% Sampling distance
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];		% Deformation fields
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
    NaN NaN NaN];

%% Downsample
% WHY DOWNSAMPLING?

matlabbatch{2}.spm.spatial.coreg.write.ref = {[templates_folder 'Template_4.nii,1']};
matlabbatch{2}.spm.spatial.coreg.write.source = {
    [output_folder 'rc1' subname '.nii,1']
    [output_folder 'rc2' subname '.nii,1']
    [output_folder 'rc3' subname '.nii,1']
    [output_folder 'rc4' subname '.nii,1']
    [output_folder 'rc5' subname '.nii,1']
    };
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'd';

%% Shoot
% Use the Shoot toolbox in SPM to register rc* to Template*
% This will output the v* (velocity) y* (exponentiated transform) j* (jacobian) files

matlabbatch{3}.spm.tools.shoot.warp1.images = {
    {[output_folder 'drc1' subname '.nii,1']}
    {[output_folder 'drc2' subname '.nii,1']}
    }';
matlabbatch{3}.spm.tools.shoot.warp1.templates = {
    [templates_folder 'Template_0.nii']
    [templates_folder 'Template_1.nii']
    [templates_folder 'Template_2.nii']
    [templates_folder 'Template_3.nii']
    [templates_folder 'Template_4.nii']
    };

spm('defaults', 'PET');
spm_jobman('run', matlabbatch);

% Delete unused files
delete([output_folder 'BiasField_' subname '.nii']);
delete([output_folder subname '_seg8.mat']);
delete([output_folder 'm' subname '.nii']);
for i1 = 2:5
    delete([output_folder 'c' num2str(i1) subname '.nii']);
    delete([output_folder 'rc' num2str(i1) subname '.nii']);
    delete([output_folder 'drc' num2str(i1) subname '.nii']);
end


%% We also need the iv* (inverse velocity)

v_dr_to_template = [output_folder 'v_drc1' subname '_Template.nii'];

nii_v_t2dr         = nifti(v_dr_to_template);
fullsize           = nii_v_t2dr.dat.dim;
nii_v_t2dr.dat.dim = nii_v_t2dr.dat.dim([1 2 3 5]);

dft  = spm_shoot_defaults;
[~,~,iv] = spm_shoot3d(single(nii_v_t2dr.dat()), [2 2 2 dft.rparam], Inf);

iv = reshape(iv, fullsize);
nii_v_t2dr.dat.dim = size(iv);

nii_v_t2dr.dat.fname = [output_folder 'iv_drc1' subname '_Template.nii'];
create(nii_v_t2dr)
nii_v_t2dr.dat(:,:,:,:,:) = iv;


%% Affine maps
% 
 mri		= [output_folder filesep subj_file];				% MRI
 mri0	= [orig_folder filesep subj_file];					% MRI0 (before alignement to MNI)
 
% crtx_nat = [output_folder 'mesh_cortex_native.gii'];
% crtx_tpl = [output_folder 'mesh_cortex_template.gii'];
% 
% % Read NIFTIs and GIIs
 nii_mri  = nifti(mri);
 nii_mri0 = nifti(mri0);
 nii_tpl  = nifti(tpl);
% gii_crtx_tpl = gifti(crtx_tpl);
% gii_crtx_nat = gifti(crtx_nat);
% 
% % Rigid mapping [mni_gareth (mm) -> native_gareth (mm)]
 tpl2native = nii_mri0.mat/nii_mri.mat;
 native2tpl = inv(tpl2native);

 save([output_folder 'affine.mat'], 'tpl2native', 'native2tpl')


% 
% % Affine mapping [mni_gareth (vox) -> dartel_gareth (vox)]
% gii_crtx_tpl2nat = transform_mesh(gii_crtx_tpl, tpl2native);
% view_mesh_volume(gii_crtx_tpl2nat, nii_mri0);	% this one...is original brought back from template space
% view_mesh_volume(gii_crtx_nat, nii_mri0);		% ...and this is the original
% view_mesh_volume(gii_crtx_tpl, nii_tpl);% new template space
% 
% figure; hold on;
% plot3(gii_crtx_tpl2nat.vertices(:,1),gii_crtx_tpl2nat.vertices(:,2),gii_crtx_tpl2nat.vertices(:,3),'ro')
% plot3(gii_crtx_nat.vertices(:,1),gii_crtx_nat.vertices(:,2),gii_crtx_nat.vertices(:,3),'g.')
% 

%% Latent code

iv_dr_to_template = [output_folder 'iv_drc1' subname '_Template.nii'];

nii_iv_t2dr         = nifti(iv_dr_to_template);
nii_iv_t2dr.dat.dim = nii_iv_t2dr.dat.dim([1 2 3 5]);

fsubspace	= [templates_folder 'subspace_scaled.nii'];
fmodel		= [templates_folder 'model_variables.mat'];
z = spm_pca_project_velocity(nii_iv_t2dr.dat(),fmodel,fsubspace);

save([output_folder 'latent_code.mat'], 'z');

%% Delete temporary files

delete([output_folder 'c1' subname '.nii']);
delete([output_folder 'rc1' subname '.nii']);
delete([output_folder 'drc1' subname '.nii']);
delete([output_folder 'iv_drc1' subname '_Template.nii']);
delete([output_folder 'j_drc1' subname '_Template.nii']);
delete([output_folder 'v_drc1' subname '_Template.nii']);
delete([output_folder 'y_drc1' subname '_Template.nii']);




