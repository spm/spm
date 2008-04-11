function vals = point_vals(obj, XYZmm, holdlist)
% returns values from all the images at points given in XYZmm
% FORMAT vals = point_vals(obj, XYZmm, holdlist)
% 
% (for the following,  I is number of images in object, N is the number
% of points to resample from)
% Input 
% obj         - object
% XYZmm       - 3xN XYZ natrix of points (in mm) 
% holdlist    - optional 1xI vector of resample hold values 
%
% Outputs
% vals        - IxN vector of values in images
% 
% $Id: point_vals.m,v 1.1 2005/04/20 15:05:36 matthewbrett Exp $ 
  
if nargin < 2
  error('Need XYZmm');
end
if nargin < 3
  holdlist = [obj.img(:).hold];
end

X=1;Y=2;Z=3;
nimgs = length(obj.img);
nvals = size(XYZmm,2);
vals = zeros(nimgs,nvals)+NaN;
if size(XYZmm,1)~=4
  XYZmm = [XYZmm(X:Z,:); ones(1,nvals)];
end
for i = 1:nimgs
  I = obj.img(i);
  XYZ = I.vol.mat\XYZmm;
  if ~mars_struct('isthere', I.vol, 'imgdata')
    vol = I.vol;
  else
    vol = I.vol.imgdata;
  end
  vals(i,:) = spm_sample_vol(vol, XYZ(X,:), XYZ(Y,:),XYZ(Z,:),[holdlist(i) ...
            I.background]);
end  
return

