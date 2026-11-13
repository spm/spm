function [mask_file, mask] = spm_mask_from_segments(mag_file, thresh, fill_holes, conn, mask_mag)

%==========================================================================
% brain extraction/masking script using spm segment tool
% it sums three probability maps from segment tool, thresholds the resulting map to obtain binary mask and fills residual holes
% Inputs:
% mag_file      - file from which to estimate the mask
% thresh        - between 0 and 1, smaller value bigger mask, threshold of the spm probability maps to create a binary map
% fill_holes    - true or false, if true it will fill small holes within the mask
% conn          - used in imfill.m, connectivity of neighbouring voxels, choose one of the following: 4, 6, 8, 18, 26
% mask_mag      - true or false, if true it will also output masked mag_file
%
% Outputs:
% mask_file     - file path from the output mask
% mask          - 3D binary mask matrix

% Barbara Dymerska
% Copyright (C) 2025 Department of Imaging Neuroscience, UCL
%==========================================================================
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

% specifying where tissue probability maps are located
TPM_file = fullfile(spm('Dir'),'tpm','TPM.nii');
TPM_V = char(TPM_file + "," + (1:6).');

% creating the structure that is required for spm_preproc8
opts.image    = spm_vol(mag_file);
opts.biasfwhm = 60;
opts.biasreg  = 0.001;
opts.tpm      = spm_load_priors8(TPM_V);
opts.lkp      = [1     2     3     3     4     4     4     5     5     5     5     6     6];
opts.reg      = [0 0.001 0.5 0.05 0.2];
opts.samp     = 3;
opts.fwhm     = 0;

Affine = spm_maff8(opts.image,3,32,opts.tpm,eye(4),'mni');
opts.Affine = spm_maff8(opts.image,3, 1,opts.tpm,Affine,'mni');

% runing the segmentation
segments =spm_preproc8(opts);

% writing out only grey, white, and CSF probability maps in native space:
write_flags = zeros(6,4);
write_flags(1:3,1) = 1;
spm_preproc_write8(segments, write_flags);


% summing up probability maps for GM, WM and CSF
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

fprintf('spm masking based on SPM segments took: %i seconds', round(toc));

