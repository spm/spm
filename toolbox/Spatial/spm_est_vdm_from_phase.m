function vdm_pha = spm_est_vdm_from_phase(vol1, vol1_mean, vol1_pha, total_readout_time, TE, PE_dir, outdir)

%==========================================================================
% Estimate Voxel Displacement Map (VDM) from phase images using ROMEO phase
% unwrapping tool.
% FORMAT:
% vdm_pha = spm_est_vdm_from_phase(vol1, vol1_mean, vol1_pha,
%                                         total_readout_time, TE, PE_dir, outdir)
%
% Input:
%   vol1                - cell array of file names of magnitude images with
%                         the same PE direction as the fMRI data to be corrected
%   vol1_mean          - file name of mean magnitude image with the same PE
%                         direction as the fMRI data to be corrected
%   vol1_pha           - cell array of file names of phase images with the
%                         same PE direction as the fMRI data to be corrected
%   total_readout_time - total readout time in seconds
%   TE                 - echo time in seconds
%   PE_dir             - phase-encoding direction (+1 or -1)
%   outdir             - output directory for temporary files
%
% Output:
%   vdm_pha            - Voxel Displacement Map (in mm) estimated from phase
%                         images
%
% Barbara Dymerska
% Copyright (C) 2025 Department of Imaging Neuroscience, UCL
%==========================================================================

if numel(vol1(:,1))~=numel(vol1_pha(:,1))
    error('You need the same nr of same PE direction phase volumes as magnitude volumes')
end

Nii_mag = zeros([size(nifti(vol1(1,:)).dat(:,:,:)) numel(vol1(:,1))]) ;
Nii_pha = zeros(size(Nii_mag)) ;

for f = 1:numel(vol1(:,1))
    Nii_mag_1vol = nifti(vol1(f,:)).dat(:,:,:);
    Nii_pha_1vol = nifti(vol1_pha(f,:)).dat(:,:,:);

    Nii_mag(:,:,:,f)= Nii_mag_1vol;
    Nii_pha(:,:,:,f)= Nii_pha_1vol;
end


% write 4D magnitude image for ROMEO phase unwrapping
%--------------------------------------------------------------------------
basename = spm_file(vol1(1,:),'basename');
oname    = char(spm_file(basename,'prefix','4D_mag_','ext','.nii'));
oname    = fullfile(outdir,oname);
Nio      = nifti;
Nio.dat  = file_array(oname,[size(Nii_mag)],'float32');
Nio.mat  = nifti(vol1(1,:)).mat;
create(Nio);
Nio.dat(:,:,:,:) = Nii_mag;

mag_file = oname;

% write 4D phase image for ROMEO phase unwrapping
%--------------------------------------------------------------------------
basename = spm_file(vol1_pha(1,:),'basename');
oname    = char(spm_file(basename,'prefix','4D_pha_','ext','.nii'));
oname    = fullfile(outdir,oname);
Nio      = nifti;
Nio.dat  = file_array(oname,[size(Nii_pha)],'float32');
Nio.mat  = nifti(vol1_pha(1,:)).mat;
create(Nio);
Nio.dat(:,:,:,:) = Nii_pha;

pha_file = oname;

% calling ROMEO phase unwrapping
ROMEO_command = sprintf('unwrapping_main(["-p", "%s", "-m", "%s", "-o", "%s", "-t", "epi", "-v", "-i", "-g", "-k", "nomask"])', pha_file, mag_file, outdir);
spm_julia('run', ROMEO_command, 'ROMEO','ArgParse', 'MriResearchTools')


pha_unwr = nifti(fullfile(outdir, 'unwrapped.nii')).dat(:,:,:,:) ;
pha_unwr = mean(pha_unwr,4) ;

vdm_pha = pha_unwr*total_readout_time/(2*pi*TE)*PE_dir ;

[~, mask] = spm_mask_from_segments(vol1_mean);

% fast VDM smoothing and extrapolation
vdm_pha = spm_smooth_extrap(vdm_pha, mask) ;


end