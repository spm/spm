function [FWHM,VRpv,R] = spm_est_smoothness(V,VM,ndf)
% Estimation of smoothness based on [residual] images
% FORMAT [FWHM,VRpv,R] = spm_est_smoothness(V,VM,[ndf]);
%
% V     - Filenames or mapped standardized residual images
% VM    - Filename of mapped mask image
% ndf   - A 2-vector, [n df], the original n & dof of the linear model
%
% FWHM  - estimated FWHM in all image directions
% VRpv  - handle of Resels per Voxel image
% R     - vector of resel counts
%__________________________________________________________________________
%
% spm_est_smoothness returns a spatial smoothness estimator based on the
% variances of the normalized spatial derivatives as described in K.
% Worsley, (1996). Inputs are a mask image and a number of standardized
% residual images, or any set of mean zero, unit variance images. Output
% is a global estimate of the smoothness expressed as the FWHM of an
% equivalent Gaussian point spread function. An estimate of resels per
% voxels (see spm_spm) is written as an image file ('RPV.img') to the
% current directory.
%
% To improve the accuracy of the smoothness estimation the error degrees
% of freedom can be supplied.  Since it is not assumed that all residual
% images are passed to this function, the full, original sample size n
% must be supplied as well.
%
% The mask image specifies voxels, used in smoothness estimation, by
% assigning them non-zero values. The dimensions, voxel sizes, orientation
% of all images must be the same. The dimensions of the images can be of
% dimensions 0, 1, 2 and 3.
%
% Note that 1-dim images (lines) must exist in the 1st dimension and
% 2-dim images (slices) in the first two dimensions. The estimated fwhm
% for any non-existing dimension is infinity.
%__________________________________________________________________________
%
% Refs:
%
% K.J. Worsley (1996). An unbiased estimator for the roughness of a
% multivariate Gaussian random field. Technical Report, Department of
% Mathematics and Statistics, McGill University
%
% S.J. Kiebel, J.B. Poline, K.J. Friston, A.P. Holmes, and K.J. Worsley.
% Robust Smoothness Estimation in Statistical Parametric Maps Using
% Standardized Residuals from the General Linear Model. NeuroImage,
% 10:756-766, 1999.
%
% S. Hayasaka, K. Phan, I. Liberzon, K.J. Worsley, T.E. Nichols (2004).
% Nonstationary cluster-size inference with random field and permutation
% methods. NeuroImage, 22:676-687, 2004.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Tom Nichols
% $Id: spm_est_smoothness.m 4136 2010-12-09 22:22:28Z guillaume $


%-Assign input arguments
%--------------------------------------------------------------------------
if nargin > 3
    spm('alert!', 'spm_est_smoothness: Wrong number of arguments');
    return;
end
if nargin < 1
    V   = spm_select(inf, '^ResI.*\.img$', 'Select residual images');
end
if nargin < 2
    VM  = spm_select(1, '^mask\.img$', 'Select mask image');
end
if nargin < 3, ndf = [NaN NaN]; end
if length(ndf) ~= 2
    error('ndf argument must be of length 2 ([n df])')
end

%-Initialise
%--------------------------------------------------------------------------
if ~isstruct(V)
    V  = spm_vol(V);
end
if ~isstruct(VM)
    VM = spm_vol(VM);
end
if any(isnan(ndf))
    ndf     = [length(V) length(V)];  % Assume full df
end
n_full = ndf(1);
edf    = ndf(2);

%-Initialise RESELS per voxel image
%--------------------------------------------------------------------------
VRpv  = struct( 'fname','RPV.img',...
    'dim',      VM.dim(1:3),...
    'dt',       [spm_type('float64') spm_platform('bigend')],...
    'mat',      VM.mat,...
    'pinfo',    [1 0 0]',...
    'descrip',  'spm_spm: resels per voxel');
VRpv   = spm_create_vol(VRpv);

%-Dimensionality of image
%--------------------------------------------------------------------------
D        = 3 - sum(VM.dim(1:3) == 1);
if D == 0
    FWHM = [Inf Inf Inf];
    R    = [0 0 0];
    return;
end

%-Find voxels within mask
%--------------------------------------------------------------------------
d          = spm_read_vols(VM);
[Ix,Iy,Iz] = ndgrid(1:VM.dim(1),1:VM.dim(2),1:VM.dim(3));
Ix = Ix(d~=0); Iy = Iy(d~=0); Iz = Iz(d~=0);

%-Compute covariance of derivatives
%--------------------------------------------------------------------------
str   = 'Spatial non-sphericity (over scans)';
fprintf('%-40s: %30s',str,'...estimating derivatives');                 %-#
spm_progress_bar('Init',100,'smoothness estimation','');

L     = zeros(size(Ix,1),D,D);
ssq   = zeros(size(Ix,1),1);
for i = 1:length(V)
    
    [d,dx,dy,dz] = spm_sample_vol(V(i),Ix,Iy,Iz,1);
    
    % sum of squares
    %----------------------------------------------------------------------
    ssq  = ssq + d.^2;
    
    % covariance of finite differences
    %----------------------------------------------------------------------
    if D >= 1
        L(:,1,1) = L(:,1,1) + dx.*dx;
    end
    if D >= 2
        L(:,1,2) = L(:,1,2) + dx.*dy;
        L(:,2,2) = L(:,2,2) + dy.*dy;
    end
    if D >= 3
        L(:,1,3) = L(:,1,3) + dx.*dz;
        L(:,2,3) = L(:,2,3) + dy.*dz;
        L(:,3,3) = L(:,3,3) + dz.*dz;
    end
    
    spm_progress_bar('Set',100*i/length(V));
    
end
spm_progress_bar('Clear')

%-Scale sum into an average (and account for DF)
% The standard result uses normalised residuals e/sqrt(RSS) and
%
%  \hat\Lambda = grad(e/sqrt(RSS))' * grad(e/sqrt(RSS))
%
% In terms of standardized residuals e/sqrt(RMS) this is
%
%  \hat\Lambda = (1/DF) * grad(e/sqrt(RMS))' * grad(e/sqrt(RMS))
%
% but both of these expressions assume that the RSS or RMS correspond to
% the full set of residuals considered.  However, spm_spm only considers
% upto MAXRES residual images.  To adapt, re-write the above as a scaled
% average over n scans
%
%  \hat\Lambda = (n/DF) * ( (1/n) * grad(e/sqrt(RMS))' * grad(e/sqrt(RMS)) )
%
% I.e. the roughness estimate \hat\Lambda is an average of outer products
% of standardized residuals (where the average is over scans), scaled by
% n/DF. Hence, we can use only a subset of scans simply by replacing this 
% last average term with an average over the subset.
%
% See Hayasaka et al, p. 678, for more on estimating roughness with
% standardized residuals (e/sqrt(RMS)) instead of normalised residuals
% (e/sqrt(RSS)). Note that the names arise from the fact that
% sqrt(RSS) = sqrt(r'*r) is norm(r), while sqrt(RMS) = sqrt(r'*r/edf)
% is the unbiased (ReML) estimate of the standard deviation.
%--------------------------------------------------------------------------
L  = L/length(V);     % Average
L  = L*(n_full/edf);  % Scale


%-Evaluate determinant (and xyz components for FWHM)
%--------------------------------------------------------------------------
if D == 1
    resel_xyz = L;
    resel_img = L;
end
if D == 2
    resel_xyz = [L(:,1,1) L(:,2,2)];
    resel_img = L(:,1,1).*L(:,2,2) - ...
                L(:,1,2).*L(:,1,2);
end
if D == 3
    resel_xyz = [L(:,1,1) L(:,2,2)  L(:,3,3)];
    resel_img = L(:,1,1).*L(:,2,2).*L(:,3,3) + ...
                L(:,1,2).*L(:,2,3).*L(:,1,3)*2 - ...
                L(:,1,1).*L(:,2,3).*L(:,2,3) - ...
                L(:,1,2).*L(:,1,2).*L(:,3,3) - ...
                L(:,1,3).*L(:,2,2).*L(:,1,3);
end    
resel_img(resel_img<0) = 0;
% Convert det(Lambda) and diag(Lambda) to units of resels
resel_img = sqrt(resel_img/(4*log(2))^D);
resel_xyz = sqrt(resel_xyz/(4*log(2)));


%-Optional mask-weighted smoothing of RPV image
%--------------------------------------------------------------------------
fwhm_vox = 0;
if any(fwhm_vox)
    if length(fwhm_vox) == 1, fwhm_vox = fwhm_vox([1 1 1]);  end
    
    % Convert resel_img at in-mask voxels to resel volume
    mask = spm_read_vols(VM) > 0;
    RPV = zeros(size(mask));
    RPV(mask) = resel_img;
    
    % Remove invalid mask voxels, i.e. edge voxels with missing derivatives
    smask = double(mask & isfinite(RPV)); % leaves mask for resel_img below
    
    % Smooth RPV volume (note that NaNs are treated as zeros in spm_smooth)
    spm_smooth(RPV, RPV, fwhm_vox);
    
    % Smooth mask and decide how far to trust smoothing-based extrapolation
    spm_smooth(smask, smask, fwhm_vox);
    infer = smask > 1e-3; % require sum of voxel's in-mask weights > 1e-3
    
    % Normalise smoothed RPV by smoothed mask
    RPV( infer) = RPV(infer) ./ smask(infer);
    RPV(~infer) = NaN; % spm_list handles remaining (unlikely) in-mask NaNs
    
    % Get data at in-mask voxels; smoothed resel_img conforms with original
    resel_img = RPV(mask);
end

%-Write Resels Per Voxel image
%--------------------------------------------------------------------------
fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...writing resels/voxel image');%-#

for i = 1:VM.dim(3)
    d = NaN(VM.dim(1:2));
    I = find(Iz == i);
    if ~isempty(I)
        d(sub2ind(VM.dim(1:2), Ix(I), Iy(I))) = resel_img(I);
    end
    VRpv = spm_write_plane(VRpv, d, i);
end


%-(unbiased) RESEL estimator and Global equivalent FWHM
% where we desire FWHM with components proportional to 1./mean(resel_xyz),
% but scaled so prod(1./FWHM) agrees with (the unbiased) mean(resel_img).
%--------------------------------------------------------------------------
i     = isnan(ssq) | ssq < sqrt(eps);
resel_img = mean(resel_img(~i,:));
resel_xyz = mean(resel_xyz(~i,:));

RESEL = resel_img^(1/D)*(resel_xyz/prod(resel_xyz)^(1/D));
FWHM  = full(sparse(1,1:D,1./RESEL,1,3));
FWHM(~FWHM) = 1;

%-resel counts for search volume (defined by mask)
%--------------------------------------------------------------------------
% R0   = spm_resels_vol(VM,[1 1 1])';
% R    = R0.*(resel.^([0:3]/3));
% OR
R      = spm_resels_vol(VM,FWHM)';


fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');               %-#
