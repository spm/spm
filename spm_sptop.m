function [K] = spm_sptop(sigma,q)
% returns a sparse Toeplitz convolution matrix
% FORMAT [K] = spm_sptop(sigma,q)
% sigma - of Gaussian kernel K (or kernel itself)
% K     - q x q sparse convolution matrix
%_______________________________________________________________________
%
% returns a q x q sparse Toeplitz Gaussian convolution matrix of sd sigma
%
%_______________________________________________________________________
% %W% Karl Friston %E%


%-----------------------------------------------------------------------
if ~sigma; K = speye(q); return; end

%-----------------------------------------------------------------------
if length(sigma) == 1

	E  = ceil(3*sigma);
	x  = [-E:E];
	k  = exp(-x.^2/(2*sigma^2));
else
	k  = sigma;
	x  = [1:length(k)] - ceil(length(k)/2);
end

K  = k(:)*ones(1,q);
j  = ones(length(k),1)*[1:q];
i  = x(:)*ones(1,q) + j;
i  = i(:);
j  = j(:);
K  = K(:);
Q  = find((i >= 1) & (i <= q));

K = sparse(i(Q),j(Q),K(Q));
