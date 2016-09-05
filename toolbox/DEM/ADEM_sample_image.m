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
% $Id: ADEM_sample_image.m 6866 2016-09-05 09:19:42Z karl $


% preliminaries
%--------------------------------------------------------------------------
global STIM

if ~isfield(STIM,'W'), STIM.W = 1/6;   end
if ~isfield(STIM,'P'), STIM.P = [0;0]; end


% retinotopic sampling
%--------------------------------------------------------------------------
dim  = size(R);
dx   = V.dim(1)/dim(1)*STIM.W;

i    = dx*((1:dim(1)) - dim(1)/2) + V.dim(1)/2  + (o(1) + STIM.P(1))*16;
j    = dx*((1:dim(2)) - dim(2)/2) + V.dim(2)/2  + (o(2) + STIM.P(2))*16;
x    = kron(ones(1,dim(2)),i);
y    = kron(j,ones(1,dim(1)));
z    = ones(1,dim(1)*dim(2));

s    = spm_sample_vol(V,x,y,z,-2);
s    = reshape(s,dim(1),dim(2)).*R;

