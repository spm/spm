function [X] = spm_conv_full(X,sx,sy)
% Hanning convolution (return full arrays)
% FORMAT [X] = spm_conv_full(X,sx,sy)
% X    - matrix
% sx   - kernel width (FWHM) in pixels
% sy   - optional non-isomorphic smoothing
%__________________________________________________________________________
%
% spm_conv_full is a one or two dimensional convolution of a matrix
% variable in working memory.  It capitalizes on the separablity of
% multidimensional convolution with a hanning kernel by using
% one-dimensional convolutions.
%__________________________________________________________________________
% Copyright (C) 1999-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_conv_full.m 7749 2019-12-05 17:05:46Z guillaume $


% assume isomorphic smoothing
%--------------------------------------------------------------------------
if nargin < 3; sy = sx; end
sx      = abs(sx);
sy      = abs(sy);
[lx,ly] = size(X);

% kernels : FWHM -> n
%--------------------------------------------------------------------------
Ex    = min([fix(sx) lx]);
kx    = spm_hanning(2*Ex + 1);
kx    = kx/sum(kx);
Ey    = min([fix(sy) ly]);
ky    = spm_hanning(2*Ey + 1);
ky    = ky/sum(ky);

% convolve
%--------------------------------------------------------------------------
if lx > 1
    for i = 1:ly
        u      = X(:,i);
        v      = [flipud(u(1:Ex)); u; flipud(u((1:Ex) + lx - Ex))];
        X(:,i) = conv(full(v),kx,'valid');
    end
end
if ly > 1
    for i = 1:lx
        u      = X(i,:);
        v      = [fliplr(u(1:Ey)) u fliplr(u((1:Ey) + ly - Ey))];
        X(i,:) = conv(full(v),ky,'valid');
    end
end
