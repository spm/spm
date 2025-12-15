function [mask_file, mask] = spm_segment_mask(mag_file, thresh, fill_holes, conn, mask_mag)

% brain extraction/masking script using spm segment tool
% it sums three probability maps from segment tool, thresholds the resulting map to obtain binary mask and fills residual holes
% Inputs:
% mag_file      - file from which to estimate the mask
% thresh        - between 0 and 1, smaller value bigger mask, threshold of the spm probability maps to create a binary map
% fill_holes    - true or false, if true it will fill small holes within the mask
% conn          - used in imfill.m, connectivity of neighbouring voxels, choose one of the following: 4, 6, 8, 18, 26
% mask_mag      - true or false, if true it will also output masked mag_file
% spm_dir       - spm directory
% output_dir    - output directory

% Barbara Dymerska

%% calling SPM segment tool
tic


if nargin < 2
    thresh = 0.5;      % masking threshold
end
if nargin < 3
    fill_holes = true; % filling holes after thresholding
end
if nargin < 4
    conn = 8;          % connectivity of neighbouring voxels
end
if nargin < 5
    mask_mag = false;  % will not mask the image
end


spm('Defaults','fMRI');
spm_jobman('initcfg');

matlabbatch{1}.spm.spatial.preproc.channel.vols = {strrep(sprintf('%s,1',mag_file),' ','')};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {sprintf('%s/tpm/TPM.nii,1', spm('Dir'))};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {sprintf('%s/tpm/TPM.nii,2', spm('Dir'))};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {sprintf('%s/tpm/TPM.nii,3', spm('Dir'))};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {sprintf('%s/tpm/TPM.nii,4', spm('Dir'))};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {sprintf('%s/tpm/TPM.nii,5', spm('Dir'))};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {sprintf('%s/tpm/TPM.nii,6', spm('Dir'))};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 0]; % saving inverse transform from MNI space to mag_file space
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];

spm_jobman('run',matlabbatch) ;
clear matlabbatch

%% summing up probability maps for GM, WM and CSF
[outdir,mag_name,~] = fileparts(mag_file) ;
mask_file = fullfile(outdir, 'mask.nii') ;

c1 = nifti(fullfile(outdir, sprintf('c1%s.nii', mag_name))) ;
c2 = nifti(fullfile(outdir, sprintf('c2%s.nii', mag_name))) ;
c3 = nifti(fullfile(outdir, sprintf('c3%s.nii', mag_name))) ;

mask(:,:,:) = c1.dat(:,:,:) + c2.dat(:,:,:) + c3.dat(:,:,:) ;
mask= imbinarize(mask,thresh) ;

if fill_holes
    mask = imfill(mask,conn,'holes') ;
else
    fprintf('not filling holes \n')
end

mask_obj = nifti(mag_file);
mask_obj.dat.fname = fullfile(outdir, 'mask.nii');
create(mask_obj)
mask_obj.dat(:,:,:) = mask ;

if mask_mag
    mag_masked_obj = nifti(mag_file) ;
    mag_masked = mag_masked_obj.dat(:,:,:).*mask(:,:,:) ;
    mag_masked_obj.dat.fname = fullfile(outdir, sprintf('%s_masked.nii',mag_name));
    create(mag_masked_obj)
    mag_masked_obj.dat(:,:,:) = mag_masked ;
end


fprintf('spm masking based on segments took: %i seconds', round(toc));

