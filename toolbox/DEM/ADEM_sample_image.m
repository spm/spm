function [s] = ADEM_sample_image(V,o,R)
% samples a (memory mapped) image at displacement o
% FORMAT [s] = ADEM_sample_image(V,o,R)
% FORMAT [s] = ADEM_sample_image(o,h)
%
% V - a structure array containing image volume information
% o - coordinates of foveal sampling:
%   o(1) - oculomotor angle
%   o(2) - oculomotor angle
% R - retinal modulation (n x n)
%
% or
%
% o - coordinates of foveal sampling
% h - vector of coefficients weighting images in STIM.H{:}
%
% s - sensory sample (n x n)
% 
% requires a global variable with the following fields:
% STIM.R = contrast modulation matrix that defines resolution
% STIM.W = width of foveal sampling of an image   (default: 1/6)
% STIM.P = image position in retinal  coordinates (default: [0;0])
% STIM.B = basis functions or receptive fields    (default: 1)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_sample_image.m 6932 2016-11-16 12:11:01Z karl $


% retinotopic predictions
%--------------------------------------------------------------------------
global STIM
if nargin < 3
    s     = 0;
    for i = 1:numel(o)
        s = s + o(i)*ADEM_sample_image(STIM.H{i},V,STIM.R);
    end
    return
end

% preliminaries
%--------------------------------------------------------------------------
if ~isfield(STIM,'R'), STIM.R = ones(64,64); end
if ~isfield(STIM,'W'), STIM.W = 1/6;         end
if ~isfield(STIM,'P'), STIM.P = [0;0];       end
if ~isfield(STIM,'B'), STIM.B = 1;           end
if ~isfield(STIM,'A'), STIM.A = 512;         end

% retinotopic sampling
%--------------------------------------------------------------------------
dim = size(R);
dx  = V.dim(1)/dim(1)*STIM.W;

i   = dx*((1:dim(1)) - dim(1)/2) + V.dim(1)/2  + (o(1) + STIM.P(1))*16;
j   = dx*((1:dim(2)) - dim(2)/2) + V.dim(2)/2  + (o(2) + STIM.P(2))*16;
x   = kron(ones(1,dim(2)),i);
y   = kron(j,ones(1,dim(1)));
z   = ones(1,dim(1)*dim(2));

x   = min(max(x,1),V.dim(1));
y   = min(max(y,1),V.dim(2));

s   = spm_sample_vol(V,x,y,z,-2); s(~s) = 1;
s   = reshape(s,dim(1),dim(2)).*R;
s   = STIM.B'*s*STIM.B;

% eccentricity attenuation
%--------------------------------------------------------------------------
s   = s*exp(-o'*o/STIM.A);

