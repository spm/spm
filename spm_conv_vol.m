
% Convolves a mapped volume with a three dimensional separable function
% FORMAT spm_conv_vol(V,Q,fx,fy,fz,offsets)
% V        -  is the memory mapped volume
% Q        -  is the output filename
% fx       -  the separable form of the function in x
% fy       -  the separable form of the function in y
% fz       -  the separable form of the function in z
% offsets  - [i j k] contains the x, y and z shifts to reposition
%             the output
% iA       -  Image array (in memory) to convolve
%_______________________________________________________________________
%
% spm_conv_vol is a compiled function (see spm_conv_vol.c).
%
% Separable means that f(x,y,z) = f(x)*f(y)*f(z) (= fx*fy*fz above)
%
% The convolution assumes zero padding in x and y with truncated smoothing 
% in z.
%
% If Q is an array with the same number of elements as the volume, the
% convolved volume will be store there instead of on disk.
%
% If iA is specified, the image dimensions are gotten from the first 3
% members of V and (if iA has the corresponding number of elements) the
% image data is read from it.
% 
% See also spm_conv.m and spm_smooth.m
%
%_______________________________________________________________________
% %W% John Ashburner, Tom Nichols %E%
