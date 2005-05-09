function [img, badvals]=pr_scaletocmap(inpimg,mn,mx,cmap,lrn)
% scales image data to colormap, returning colormap indices
% FORMAT [img, badvals]=pr_scaletocmap(inpimg,mn,mx,cmap,lrn)
% 
% Inputs
% inpimg     - matrix containing image to scale
% mn         - image value that maps to first value of colormap
% mx         - image value that maps to last value of colormap
% cmap       - 3xN colormap
% lrn        - 1x3 vector, giving colormap indices that should fill:
%              - lrn(1) (L=Left) - values less than mn
%              - lrn(2) (R=Right) - values greater than mx
%              - lrn(3) (N=NaN) - NaN values
%             If lrn value is 0, then colormap values are set to 1, and
%             indices to these values are returned in badvals (below)
% 
% Output
% img        - inpimg scaled between 1 and (size(cmap, 1))
% badvals    - indices into inpimg containing values out of range, as
%              specified by lrn vector above
% 
% $Id: pr_scaletocmap.m,v 1.1 2005/04/20 15:05:00 matthewbrett Exp $

if nargin <  4
  error('Need inpimg, mn, mx, and cmap');
end

cml = size(cmap,1);

if nargin < 5
  lrn = [1 cml 0];
end

img = (inpimg-mn)/(mx-mn);  % img normalized to mn=0,mx=1
if cml==1 % values between 0 and 1 -> 1
  img(img>=0 & img<=1)=1;
else
  img = img*(cml-1)+1;
end
outvals = {img<1, img>cml, isnan(img)};
img= round(img);
badvals = zeros(size(img));
for i = 1:length(lrn)
  if lrn(i)
    img(outvals{i}) = lrn(i);
  else
    badvals = badvals | outvals{i};
    img(outvals{i}) = 1;
  end    
end
return
