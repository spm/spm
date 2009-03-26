function [fwhm,VRpv] = spm_est_smoothness(varargin)
% Estimation of smoothness based on [residual] images
% FORMAT [fwhm,VRpv] = spm_est_smoothness(VResI,[VM,ndf]);
%
% VResI - Filenames or mapped standardized residual images
% VM    - Filename of mapped mask image
% ndf   - A 2-vector, [n df], the original n & dof of the linear model
%
% fwhm  - estimated FWHM in all image directions
% VRpv  - handle of Resels per Voxel image
%_______________________________________________________________________
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
%_______________________________________________________________________
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
% S. Hayasaka, K. Phan, I. Liberzon, K.J. Worsley, T.E .Nichols (2004).
% Nonstationary cluster-size inference with random field and permutation
% methods. NeuroImage, 22:676-687, 2004.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Tom Nichols
% $Id: spm_est_smoothness.m 2957 2009-03-26 11:18:48Z guillaume $


%-Assign input arguments
%-----------------------------------------------------------------------
if nargin > 0, V   = varargin{1}; end
if nargin > 1, VM  = varargin{2}; end
if nargin > 2, ndf = varargin{3}; end
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
%-----------------------------------------------------------------------
if ~isstruct(V)
    V     = spm_vol(V);
end
if ~isstruct(VM)
    VM    = spm_vol(VM);
end
if any(isnan(ndf))
  ndf     = [length(V) length(V)];  % Assume full df
end
n_full    = ndf(1);
edf       = ndf(2);

%-Intialise RESELS per voxel image
%-----------------------------------------------------------------------
VRpv  = struct('fname','RPV.img',...
            'dim',      VM.dim(1:3),...
            'dt',       [spm_type('float64') spm_platform('bigend')],...
            'mat',      VM.mat,...
            'pinfo',    [1 0 0]',...
            'descrip',  'spm_spm: resels per voxel');
VRpv  = spm_create_vol(VRpv);


%-Dimensionality of image
%-----------------------------------------------------------------------
N     = 3 - sum(VM.dim(1:3) == 1);
if N == 0
    fwhm = [Inf Inf Inf];
    return
end

%-Find voxels within mask
%-----------------------------------------------------------------------
[x,y] = ndgrid(1:VM.dim(1), 1:VM.dim(2));
I     = []; Ix = []; Iy = []; Iz = [];
for i = 1:VM.dim(3)
    z  = i*ones(size(x));
    d  = spm_sample_vol(VM, x, y, z, 0);
    I  = find(d);
    Ix = [Ix; x(I)];
    Iy = [Iy; y(I)];
    Iz = [Iz; z(I)];
end

%-Compute variance of derivatives in all directions
%-----------------------------------------------------------------------
str   = 'Spatial non-sphericity (over scans)';
fprintf('%-40s: %30s',str,'...estimating derivatives');              %-#
spm_progress_bar('Init',100,'smoothness estimation','');

v     = zeros(size(Ix,1),N);
ssq   = zeros(size(Ix,1),1);
for i = 1:length(V) % for all residual images
    
    [d, dx, dy, dz] = spm_sample_vol(V(i), Ix, Iy, Iz, 1);
    
    if N >= 1. v(:, 1) = v(:, 1) + dx.^2;  end
    if N >= 2. v(:, 2) = v(:, 2) + dy.^2;  end
    if N >= 3, v(:, 3) = v(:, 3) + dz.^2;  end

    ssq  = ssq + d.^2;

    spm_progress_bar('Set',100*i/length(V));

end
spm_progress_bar('Clear')

%-Scale sum into an average (and account for DF)
%
% The standard result uses normalised residuals e/sqrt(RSS) and
%
%    \hat\Lambda = grad(e/sqrt(RSS))' * grad(e/sqrt(RSS))
%
% Note this function (for now) neglects the off diagonals of \Lambda.
%
% In terms of standardized residuals e/sqrt(RMS) this is
%
%    \hat\Lambda = (1/DF) * grad(e/sqrt(RMS))' * grad(e/sqrt(RMS))
%
% but both of these expressions assume that the RSS or RMS correspond to
% the full set of residuals considered.  However, spm_spm only considers
% upto MAXRES residual images.  To adapt, re-write the above as a scaled
% average over n scans
%
%    \hat\Lambda = (n/DF) * ( (1/n) * grad(e/sqrt(RMS))' * grad(e/sqrt(RMS)) )
%
% I.e. the roughness estimate \hat\Lambda is an average of outer products
% of standardized residuals (where the average is over scans), scaled by
% n/DF.
%
% Hence, we can use only a subset of scans simply by replacing this last
% average term with an average over the subset.
% 
% See Hayasaka et al, p. 678, for more on estimating roughness with
% standardized residuals (e/sqrt(RMS)) instead of normalised residuals
% (e/sqrt(RSS)).
%-----------------------------------------------------------------------
v = v / length(V); % Average 

v = v * (n_full/edf); % Scale


%-Eliminate NaN / zero variance voxels
%-----------------------------------------------------------------------
I      = find(any(isnan(v),2) | ssq<sqrt(eps));
v(I,:) = []; Ix(I) = []; Iy(I) = []; Iz(I) = [];


%-Write Resels Per Voxel image
%-----------------------------------------------------------------------
fprintf('\r%-40s: %30s\n',str,'...writing resels/voxel image');      %-#

resel_img_3 = sqrt(v./(4*log(2)));
resel_img   = prod(resel_img_3,2);
for  i = 1:VM.dim(3)
    d  = NaN(VM.dim(1:2));
    I  = find(Iz == i);
    if ~isempty(I)
        d(sub2ind(VM.dim(1:2), Ix(I), Iy(I))) = resel_img(I);
    end
    VRpv = spm_write_plane(VRpv, d, i);
end

%-Global equivalent FWHM
%-----------------------------------------------------------------------
resel    = mean(resel_img);
resel_3  = mean(resel_img_3);
FWHM     = 1 ./ resel_3;

%-Bias correction {prod(1/FWHM) = (unbiased) RESEL estimator}
%-----------------------------------------------------------------------
% Variable 'resel' is unbiased and is actual quantity that determines
% RFT P-values (see spm_resels_vol and use of SPM.xVol.R), but only FWHM
% gets returned from this function. Hence the following bias correction
% ensures that prod(1./FWHM) equals resel. 
FWHM     = FWHM / prod(FWHM)^(1/N) * (1/resel).^(1/N);

%-Carefully fill-in accounting for dimension
fwhm     = [Inf Inf Inf];
fwhm(1:length(FWHM)) = FWHM;
