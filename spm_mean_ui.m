function spm_mean
% promts for a series of images and averages them
% FORMAT spm_mean
%_______________________________________________________________________
%
% spm_mean simply averages a set of images to produce a mean image
% that is written to "mean.img" (in the current directory).
% Any spatial differences described by the `mat' files are ignored, and
% the resultant mean is given the spatial orientation of the first image.
%_______________________________________________________________________
% %W% John Ashburner %E%

%-Select images
%-----------------------------------------------------------------------
V=spm_vol(spm_get(Inf,'.img','Select images to be averaged'));

%-Compute mean and write headers etc.
%-----------------------------------------------------------------------
VO         = V(1);
VO.fname   = 'mean.img';
VO.descrip = 'Mean image';
VO.pinfo   = [1.0 0 0]';
VO.dim(4)  = 4;
spm_create_image(VO);
for i=1:prod(size(V)), V(i).pinfo(1:2,:) = V(i).pinfo(1:2,:)/prod(size(V)); end;
VO.pinfo(1,1) = spm_add(V,VO);
spm_create_image(VO);
