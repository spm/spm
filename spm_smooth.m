function spm_smooth(P,Q,s,Vi)
% 3 dimensional convolution of an image
% FORMAT spm_smooth(P,Q,S,Vi)
% P  - image to be smoothed
% Q  - filename for smoothed image
% S  - [sx sy sz] Guassian filter width {FWHM} in mm
% Vi - Volume info (only used if P is an array)
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
%__________________________________________________________________________
% %W% Anon, Tom Nichols %E%

%----------------------------------------------------------------------------
if length(s) == 1; s = [s s s]; end


% read and write header if we're working with files
%----------------------------------------------------------------------------
if isstr(P)
    [DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P);
else
    if nargin < 4
	error(['spm_smooth: Must specify image descriptor vector ', ...
	    'if smoothing image from RAM']); end
    DIM = Vi(1:3)';
    VOX = Vi(4:6)';
    SCALE = Vi(7);
    TYPE = Vi(8);
    OFFSET = Vi(9);
    ORIGIN = Vi(10:12);
    DESCRIP = '';
end

if isstr(Q)
    spm_hwrite(Q,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,[DESCRIP ' -conv']); end

    
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

if isstr(P)
    V  = spm_map(P);
    spm_conv_vol(V,Q,x,y,z,(-[i,j,k]));
    spm_unmap(V);
else
    spm_conv_vol(Vi,Q,x,y,z,(-[i,j,k]),P);
end




