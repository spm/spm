function [O] = spm_conv(X,sx,sy)
% Gaussian convolution
% FORMAT [O] = spm_conv(X,sx,[sy]);
% X    - matrix
% sx   - kernel width (FWHM) in pixels
% sy   - optional non-isomorphic smoothing
%___________________________________________________________________________
%
% spm_conv is a one or two dimensional convolution of a matrix variable
% in working memory.  It capitalizes on the sparisity structure of the
% problem and the separablity of multidimensional convolution with a Gaussian
% kernel by using one-dimensional convolutions and kernels that are
% restricted to non near-zero values
%
%__________________________________________________________________________
% %W% %E%

% assume isomorphic smoothing
%---------------------------------------------------------------------------
if nargin < 3; sy = sx; end

% FWHM -> sigma
%---------------------------------------------------------------------------
sx    = sx/sqrt(8*log(2));
sy    = sy/sqrt(8*log(2));

% create and mulitply by convolution matrices
%---------------------------------------------------------------------------
[lx ly] = size(X);
K       = spm_sptop(sx,lx);
O       = K*X;
if (ly ~= lx) | (sy ~= sx)
	K  = spm_sptop(sy,ly);
end
O       = full(K*O')';
