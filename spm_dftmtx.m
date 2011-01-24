function C = spm_dftmtx(N,K,a)
% Creates basis functions for Discrete Cosine Transform.
% FORMAT C = spm_dftmtx(N,K,a)
%
% N - dimension
% K - order
% a - number of (1/2)5Hz frequency steps (default a = 2)
%__________________________________________________________________________
% spm_dftmtx creates a matrix for the first few basis functions of a one
% dimensional discrete Fourier transform.
%
% See:    Fundamentals of Digital Image Processing (p 150-154).
%         Anil K. Jain 1989.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dftmtx.m 4170 2011-01-24 18:37:42Z karl $
 
 
% initialise
%--------------------------------------------------------------------------
try, a; catch, a = 2; end
 
C     = ones(N,1);
x     = (0:(N - 1))/N;
for k = 1:K
    C(:,end + 1) = sin(a*pi*k*x);
    C(:,end + 1) = cos(a*pi*k*x);
end
