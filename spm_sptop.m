function [K] = spm_sptop(sigma,q)
% returns a sparse Toeplitz convolution matrix
% FORMAT [K] = spm_sptop(sigma,q)
% sigma - of Gaussian kernel K (or kernel itself)
% K     - q x q sparse convolution matrix
%_______________________________________________________________________
%
% returns a q x q sparse convolution matrix
%
%_______________________________________________________________________
% %W% Karl Friston %E%

% if simgma = 0 return identity matrix, if q = 1 return 1.
%-----------------------------------------------------------------------
if ~any(sigma); K = speye(q); return; end
if q == 1;      K = 1;        return; end

% otherwise get kernel function
%-----------------------------------------------------------------------
if length(sigma) == 1

	E  = ceil(3*sigma);
	x  = [-E:E];
	k  = exp(-x.^2/(2*sigma^2));
else
	k  = sigma;
	x  = [1:length(k)] - 1;
end

% and create convolution matrix
%-----------------------------------------------------------------------
K  = k(:)*ones(1,q);
j  = ones(length(k),1)*[1:q];
i  = x(:)*ones(1,q) + j;

% setting the row-wise sum to unity
%-----------------------------------------------------------------------
i  = i(:);
j  = j(:);
K  = K(:);
Q  = find((i >= 1) & (i <= q));
K  = sparse(i(Q),j(Q),K(Q));
Q  = sum(K');
K  = inv(diag(Q + (~Q)))*K;
