function [s]= ADEM_sample_image(V,o,R)
% samples a (memory mapped) image at displacement o
% FORMAT [s]= ADEM_sample_image(V,o,R)
%
% V - a structure array containing image volume information
% o - hidden states:
%   o(1) - oculomotor angle
%   o(2) - oculomotor angle
% R - retinal modulation (n x n)
%
% s - sensory input (n x n)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_sample_image.m 4595 2011-12-19 13:06:22Z karl $
 

% retinotopic sampling
%--------------------------------------------------------------------------
dim  = size(R);
dx   = V.dim(1)/dim(1)/6;

i    = dx*((1:dim(1)) - dim(1)/2) + V.dim(1)/2  + o(1)*16;
j    = dx*((1:dim(2)) - dim(2)/2) + V.dim(2)/2  + o(2)*16;
x    = kron(ones(1,dim(2)),i);
y    = kron(j,ones(1,dim(1)));
z    = ones(1,dim(1)*dim(2));

s    = spm_sample_vol(V,x,y,z,-2);
s    = reshape(s,dim(1),dim(2)).*R;

