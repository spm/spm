function pm = pm_smooth_phasemap(pm,angvar,vxs,fwhm)
% Performs a weighted (by 1/angvar) gaussian smoothing
% of a phasemap.
% FORMAT: pm = pm_smooth_phasemap(pm,angvar,vxs,fwhm)
%
% Input:
% pm         : Phase-map
% angvar     : Map of uncertainty of the angular estimate.
% vxs        : Voxel sizes (mm) in the three directions.
% fwhm       : FWHM (mm) of gaussian kernel for the three 
%              directions (or scalar for isotropic kernel).
%__________________________________________________________
% Jesper Andersson 16/10-03
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: pm_smooth_phasemap.m 6791 2016-04-28 14:47:20Z john $

if nargin ~= 4 || nargout ~= 1
   help pm_smooth_phasemap
   pm = [];
   return;
end
if length(fwhm) == 1
   fwhm = fwhm * [1 1 1];
end
 
%
% Determine kernel. The code has been knicked
% (and changed a bit) from spm_smooth.
%
fwhm  = fwhm./vxs;                  % voxel anisotropy
fwhm  = max(fwhm,ones(size(fwhm))); % lower bound on FWHM
s     = fwhm/sqrt(8*log(2));        % FWHM -> Gaussian parameter

x  = round(3*s(1)); x = [-x:x];
y  = round(3*s(2)); y = [-y:y];
z  = round(3*s(3)); z = [-z:z];
x  = exp(-(x).^2/(2*(s(1)).^2)); 
y  = exp(-(y).^2/(2*(s(2)).^2)); 
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;

q1 = zeros(size(pm));
q2 = zeros(size(pm));
spm_conv_vol(pm./angvar,q1,x,y,z,-[i,j,k]);
spm_conv_vol( 1./angvar,q2,x,y,z,-[i,j,k]);
pm = q1./q2;

return

