
% convolves a mapped volume with a three dimensional separable function
% FORMAT spm_conv_vol(V,Q,fx,fy,fz,offsets)
% V    -  is the memory mapped volume
% Q    -  is the output filename
% fx   -  the separable form of the function in x
% fy   -  the separable form of the function in y
% fz   -  the separable form of the function in z
% offsets  - [i j k] contains the x, y and z shifts to reposition the output
%____________________________________________________________________________
%
% spm_conv_vol is a compiled function (see spm_conv_vol.c).  Separable means
% that f(x,y,z) = f(x)*f(y)*f(z) (= fx*fy*fz above)
% The convolution assumes zero padding in x and y with truncated smoothing 
% in z.
%
% see also spm_conv.m and spm_smooth.m
%
%__________________________________________________________________________
% %W% %E%
