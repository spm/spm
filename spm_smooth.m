function spm_smooth(P,Q,s)
% 3 dimensional convolution of an image
% FORMAT spm_smooth(P,Q,S)
% P  - image to be smoothed
% Q  - filename for smoothed image
% S  - [sx sy sz] Guassian filter width {FWHM} in mm
%____________________________________________________________________________
%
% spm_smooth is used to smooth or convolve images in a file.
%
% The sum of kernel coeficients are set to unity.  Boundary
% conditions assume data does not exist outside the image in z (i.e.
% the kernel is truncated in z at the boundaries of the image space). S
% can be a vector of 3 FWHM values that specifiy an anisotropic
% smoothing.  If S is a scalar isotropic smoothing is implemented.
%
%__________________________________________________________________________
% %W% %E%

%----------------------------------------------------------------------------
if length(s) == 1; s = [s s s]; end

% write header
%----------------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P);
spm_hwrite(Q,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,[DESCRIP ' -conv']);

% compute parameters for spm_conv_vol
%----------------------------------------------------------------------------
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

V  = spm_map(P);
spm_conv_vol(V,Q,x,y,z,(-[i,j,k]));
spm_unmap_vol(V);



