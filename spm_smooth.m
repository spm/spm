function spm_smooth(P,Q,s)
% 3 dimensional convolution of an image
% FORMAT spm_smooth(P,Q,S,Vi)
% P  - image to be smoothed
% Q  - filename for smoothed image
% S  - [sx sy sz] Guassian filter width {FWHM} in mm
%____________________________________________________________________________
%
% spm_smooth is used to smooth or convolve images in a file (maybe).
%
% The sum of kernel coeficients are set to unity.  Boundary
% conditions assume data does not exist outside the image in z (i.e.
% the kernel is truncated in z at the boundaries of the image space). S
% can be a vector of 3 FWHM values that specifiy an anisotropic
% smoothing.  If S is a scalar isotropic smoothing is implemented.
%
% If P is not a string, it is taken to be volume data whose dimensions
% and voxel sizes are specified in Vi.
%
% If Q is not a string, it is used as the destination of the smoothed
% image.  It must already be defined with the same number of elements 
% as the image. 
%
%_______________________________________________________________________
% %W% John Ashburner, Tom Nichols %E%

%-----------------------------------------------------------------------
if length(s) == 1; s = [s s s]; end


% read and write header if we're working with files
%-----------------------------------------------------------------------
if isstr(P)
	V   = spm_vol(P);
	VOX = sqrt(sum(V.mat(1:3,1:3).^2));
else
    if nargin < 4
	error(['spm_smooth: Must specify image descriptor vector ', ...
	    'if smoothing image from RAM']); end
    VOX = [1 1 1];
end

if isstr(Q) & isstr(P),
	VO         = V;
	VO.fname   = Q;
	VO.descrip = sprintf('%s -conv (%g,%g,%g)',DESCRIP,s);
	spm_create_image(VO);
end


% compute parameters for spm_conv_vol
%-----------------------------------------------------------------------
s  = s./VOX;					% voxel anisotropy
s  = max(s,ones(size(s)));			% lower bound on FWHM
s  = s/sqrt(8*log(2));				% FWHM -> Gaussian parameter

x  = round(6*s(1)); x = [-x:x];
y  = round(6*s(2)); y = [-y:y];
z  = round(6*s(3)); z = [-z:z];
x  = exp(-(x).^2/(2*(s(1)).^2)); 
y  = exp(-(y).^2/(2*(s(2)).^2)); 
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;

if isstr(P)
    spm_conv_vol(spm_vol(P),Q,x,y,z,-[i,j,k]);
else
    spm_conv_vol(P,Q,x,y,z,-[i,j,k]);
end
