function [K] = spm_sptop(sigma,q,c)
% Sparse Toeplitz convolution matrix given convolution kernel
% FORMAT [K] = spm_sptop(sigma,q,c)
%
% sigma - of Gaussian kernel K (or kernel itself)
% q     - order of matrix
% c     - kernel index at t = 0 {default c = length(sigma)/2) 
% K     - q x q sparse convolution matrix
%__________________________________________________________________________
%
% Returns a q x q sparse convolution matrix.  If sigma is a scalar then
% a symmetrical Gaussian convolution matrix is returned with kernel width
% = sigma.  If sigma is a vector than sigma constitutes the kernel.  To
% obtain an asymmetrical convolution matrix (i.e. implement a phase shift
% set c = 1.
%
% Boundary handling: The row-wise sum of K is set to unity (kernel truncation)
%
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1996-2022 Wellcome Centre for Human Neuroimaging


% if sigma = 0, return identity matrix; if q = 1, return 1.
%--------------------------------------------------------------------------
if ~any(sigma); K = speye(q); return; end
if q == 1;      K = 1;        return; end

% otherwise get kernel function
%--------------------------------------------------------------------------
if length(sigma) == 1
    E  = ceil(3*sigma);
    x  = [-E:E];
    k  = exp(-x.^2/(2*sigma^2));
else
    if nargin < 3
        c = length(sigma)/2;
    end
    E  = length(sigma);
    x  = [1:E] - fix(c);
    k  = sigma;
end

% and create convolution matrix
%--------------------------------------------------------------------------
K  = k(:)*ones(1,q);
j  = ones(length(k),1)*[1:q];
i  = x(:)*ones(1,q) + j;

% setting the row-wise sum to unity
%--------------------------------------------------------------------------
Q  = find((i >= 1) & (i <= q));
K  = sparse(i(Q),j(Q),K(Q));
Q  = sum(K');
K  = inv(diag(Q + (~Q)))*K;
