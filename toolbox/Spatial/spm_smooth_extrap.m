function image = spm_smooth_extrap(image, mask, k, Niter)
% Fast convolution Gaussian smoothing in 3D plus extrapolation to missing data
% FORMAT v = spm_smooth_extrap(image, mask, k, Niter)
% Inputs:
%   image: 3D image to be smoothed  
%   mask: 3D binary mask defining voxels with valid data
%   k: Scalar smoothing kernel in image domain specified in voxels
%   Niter: number of iteration steps used for extrapolation
% Output:
%   image: Smoothed 3D image
%
%   Barbara Dymerska, Oliver Josephs
%   Copyright (C) 2025 Department of Imaging Neuroscience, UCL
%==========================================================================

eps = 1e-5;

if nargin <2
    mask = ones(size(image)) ;
end

if nargin <3
    k = 2 ;
end

if nargin <4
    Niter = 10 ;
end


image_orig = image ;
mask_orig = mask ;

for iter = 1:Niter
image = spm_smooth_fft(image.*mask, k) ;
% slightly increased smoothing kernel to get smooth edges of extrapolated region:
mask = spm_smooth_fft(single(mask), k+0.1); 

image = image ./ max(mask, eps);
if iter<Niter
image(mask_orig) = image_orig(mask_orig) ;
end

end

end

%% ===========================================================================
function v = spm_smooth_fft(v, k)

dims=size(v);

% Fourier Transform in all 3 dims:
v=fftn(v) ;

% Coordinate grids
x = (1:dims(1)) - dims(1)/2 - 1;
y = (1:dims(2)) - dims(2)/2 - 1;
z = (1:dims(3)) - dims(3)/2 - 1;

% Form a real-space Gaussian and FFT to k-space
ga=exp(ifftshift(single(-( x'.^2 + y.^2 + shiftdim(z.^2, -1) ) / k^2)));
ga=fftn(ga);

% Unit kernel volume in image space <=> unit amplitude at "dc point" in k space
ga=ga/ga(1);

% Apply smoothing kernel via efficient multiplication in k-space:
v=v.*ga;

% Return to image space:
v = real(ifftn(v)) ;

end