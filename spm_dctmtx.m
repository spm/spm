% Creates basis functions for Discrete Cosine Transform.
% FORMAT C = spm_dctmtx(N,K,n);
% N - dimension
% K - order
% n - optional points to sample
%____________________________________________________________________________
% spm_dctmtx creates a matrix for the first few basis functions of a one
% dimensional discrete cosine transform.
% See:    Fundamentals of Digital Image Processing (p 150-154).
%         Anil K. Jain 1989.
%

% %W% John Ashburner MRCCU/FIL %E%

function C = spm_dctmtx(N,K,n)

if (nargin == 2)
	n = (0:(N-1))';
else
	n = n(:);
end

C = zeros(size(n,1),K);
C(:,1)=ones(size(n,1),1)/sqrt(N);

for k=2:K
	C(:,k) = sqrt(2/N)*cos(pi*(2*n+1)*(k-1)/(2*N));
end


