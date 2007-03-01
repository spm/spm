function mask = pm_mask(angvar,mthres,ndil)
%
% Creating a mask that will determine how far
% we should proceed with phase unwrapping.
% FORMAT: mask = pm_mask(angvar,mthrea,ndil)
%
% Input:
% angvar     : Map of variance of angle estimate.
% mthres     : Threshold for variance beyond which
%              phase unwrapping is considered too
%              uncertain. Default value (pi^2)/6
%              is half the variance of a U[-pi,pi]
%              distribution.
% ndil       : We can optionally specify a no. of 
%              erodes-dilates to apply to the mask
%              in order to exclude areas connected
%              only by thin bridges to the rest of
%              the brain.
%
% Output:
% mask       : Well...
%
%_____________________________________________________
% Jesper Andersson 26/9-03.

if nargin < 2
   mthres = (pi^2)/6;
end
if nargin < 3
   ndil = 0;
end

%
% Threshold map of angular variance and pick
% largest connected component.
%
mask = double(angvar < mthres);
[lmask,num] = ip_bwlabel(mask,6);
n = histc(lmask(:),[0:num]+0.5);
[mv,mi] = max(n);
indx = lmask(:)==mi;
mask(~indx) = 0;

if ndil
   dmask = mask;
   for i=1:ndil
      dmask = ip_erode(dmask);
   end
   [lmask,num] = ip_bwlabel(dmask,6);
   n = histc(lmask(:),[0:num]+0.5);
   [mv,mi] = max(n);
   indx = lmask(:)==mi;
   dmask(~indx) = 0;
   for i=1:ndil
      dmask = ip_dilate(dmask);
   end
   mask = mask.*dmask;
end

return
 
