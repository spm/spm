function [X] = spm_conv(X,sx,sy)
% Gaussian convolution
% FORMAT [X] = spm_conv(X,sx[,sy])
% X    - matrix
% sx   - kernel width (FWHM) in pixels
% sy   - optional non-isomorphic smoothing
%__________________________________________________________________________
%
% spm_conv is a one or two dimensional convolution of a matrix variable in
% working memory.  It capitalizes on the sparsity structure of the problem
% and the separablity of multidimensional convolution with a Gaussian
% kernel by using one-dimensional convolutions and kernels that are
% restricted to non near-zero values.
%__________________________________________________________________________
% Copyright (C) 1999-2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_conv.m 7714 2019-11-26 11:25:50Z spm $


% assume isomorphic smoothing
%--------------------------------------------------------------------------
if nargin < 3; sy = sx; end
sx      = abs(sx);
sy      = abs(sy);
[lx,ly] = size(X);

% FWHM -> sigma
%--------------------------------------------------------------------------
sx    = sx/sqrt(8*log(2)) + eps;
sy    = sy/sqrt(8*log(2)) + eps;

% kernels
%--------------------------------------------------------------------------
Ex    = min([fix(6*sx) lx]);
x     = -Ex:Ex;
kx    = exp(-x.^2/(2*sx^2));
kx    = kx/sum(kx);
Ey    = min([fix(6*sy) ly]);
y     = -Ey:Ey;
ky    = exp(-y.^2/(2*sy^2));
ky    = ky/sum(ky);

% convolve
%--------------------------------------------------------------------------
if spm_check_version('matlab','7.5') > 0
    shape = {'valid'};
else
    shape = {};
end
if lx > 1
    for i = 1:ly
        u      = X(:,i);
        v      = [flipud(u(1:Ex)); u; flipud(u((1:Ex) + lx - Ex))];
        if isempty(shape)
            V  = sparse(conv(full(v),kx));
            X(:,i) = V([1:lx] + 2*Ex);
        else
            X(:,i) = sparse(conv(full(v),kx,shape{:}));
        end
    end
end
if ly > 1
    for i = 1:lx
        u      = X(i,:);
        v      = [fliplr(u(1:Ey)) u fliplr(u((1:Ey) + ly - Ey))];
        if isempty(shape)
            V  = sparse(conv(full(v),ky));
            X(i,:) = V([1:ly] + 2*Ey);
        else
            X(i,:) = sparse(conv(full(v),ky,shape{:}));
        end
    end
end
