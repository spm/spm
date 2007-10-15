function spm_smooth(P,Q,s,dtype)
% 3 dimensional convolution of an image
% FORMAT spm_smooth(P,Q,S,dtype)
% P  - image to be smoothed
% Q  - filename for smoothed image
% S  - [sx sy sz] Guassian filter width {FWHM} in mm
% dtype - datatype
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
% If Q is not a string, it is used as the destination of the smoothed
% image.  It must already be defined with the same number of elements
% as the image.
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner & Tom Nichols
% $Id: spm_smooth.m 946 2007-10-15 16:36:06Z john $


%-----------------------------------------------------------------------
if length(s) == 1; s = [s s s]; end
if nargin<4, dtype = 0; end;

if ischar(P),

    % Peter Stiers' bug fix for 4D data
    if size(P,1)>1,
        for j=1:size(P,1),
            spm_smooth(deblank(P(j,:)),Q,s,dtype);
        end
        return;
    end

    P = spm_vol(P);
end;
if isstruct(P),
    VOX = sqrt(sum(P.mat(1:3,1:3).^2));
else
    VOX = [1 1 1];
end;

if ischar(Q) && isstruct(P),
    [pth,nam,ext,num] = spm_fileparts(Q);
    q         = fullfile(pth,[nam,ext]);
    Q         = P;
    Q.fname   = q;
    if ~isempty(num),
        Q.n       = str2num(num);
    end;
    if ~isfield(Q,'descrip'), Q.descrip = sprintf('SPM compatible'); end;
    Q.descrip = sprintf('%s - conv(%g,%g,%g)',Q.descrip, s);

    if dtype~=0, % Need to figure out some rescaling.
        Q.dt(1) = dtype;
        if ~isfinite(spm_type(Q.dt(1),'maxval')),
            Q.pinfo = [1 0 0]'; % float or double, so scalefactor of 1
        else
            % Need to determine the range of intensities
            if isfinite(spm_type(P.dt(1),'maxval')),
                % Integer types have a defined maximum value
                maxv = spm_type(P.dt(1),'maxval')*P.pinfo(1) + P.pinfo(2);
            else
                % Need to find the max and min values in original image
                mx = -Inf;
                mn =  Inf;
                for pl=1:P.dim(3),
                    tmp = spm_slice_vol(P,spm_matrix([0 0 pl]),P.dim(1:2),0);
                    tmp = tmp(isfinite(tmp));
                    mx  = max(max(tmp),mx);
                    mn  = min(min(tmp),mn);
                end
                maxv = max(mx,-mn);
            end
            sf      = maxv/spm_type(Q.dt(1),'maxval');
            Q.pinfo = [sf 0 0]';
        end
    end
end

% compute parameters for spm_conv_vol
%-----------------------------------------------------------------------
s  = s./VOX;				     	% voxel anisotropy
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


if isstruct(Q), Q = spm_create_vol(Q); end;
spm_conv_vol(P,Q,x,y,z,-[i,j,k]);
