function ima = spm_dilate(varargin)
% Perform a dilation on an image (2D or 3D) 
% It uses either the supplied kernel or a standard 6-connectivity kernel.
% FORMAT: ima = spm_dilate(ima)
% or
% FORMAT: ima = spm_dilate(ima,kernel)
%
% Input:
% ima    : 2 or 3D image
% kernel : (Optional) voxel values in ima are replaced by the 
%          maximum value in a neighbourhood defined by kernel.
%          The "standard" dilation operation (in 2D) is realised
%          using the kernel
%          0 1 0
%          1 1 1
%          0 1 0
%
% Output:
% ima    : Dilated image.
%
% The functionality of this routine has been modelled on the function
% imdilate from the MATLAB Image processing toolbox. It doesn't (yet)
% have a support function such as strel to help the user to define
% kernels (you have to do it yourself if you want anything above
% 6-connectivty) and it doesnt do the clever structuring element
% decomposition that strel does (and imdilate uses). That should
% in principle mean that spm_dilate is slower than imdilate, but
% at least for small (typical) kernels it is actually more than
% twice as fast.
% The actual job is done by spm_dilate_erode.c that serves both
% spm_dilate.m and spm_erode.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id: spm_dilate.m 4310 2011-04-18 16:07:35Z guillaume $


if exist('spm_dilate_erode','file')~=3 
   error('spm_dilate_erode.c not compiled - see Makefile');
end
   
if nargin > 1
   kernel = varargin{2};
else
   if length(size(varargin{1})) == 2
      kernel = [0 1 0; 1 1 1; 0 1 0];
   elseif length(size(varargin{1})) == 3
      kernel = cat(3,[0 0 0; 0 1 0; 0 0 0],[0 1 0; 1 1 1; 0 1 0],[0 0 0; 0 1 0; 0 0 0]);
   else
      error('Input ima must be 2- or 3-dimensional');
   end
end

ima = spm_dilate_erode(varargin{1},kernel,'dilate');
