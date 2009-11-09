function [x] = spm_samp_gauss (m, C, N, dC, vC)
% Sample from a Gaussian PDF
% FORMAT [x] = spm_samp_gauss (m, C, N, dC, vC)
% m     [d x 1] mean
% C     [d x d] covar
% N     Number of samples
% dC    diagonalised C [d x 1]
% vC    eigenvectors of C [d x d]
%
% x     [N x d] matrix of samples
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_samp_gauss.m 3543 2009-11-09 09:40:46Z maria $

d = size(C, 1);
m = reshape(m, 1, d);   % Ensure that m is a row vector

if nargin > 3
    deig = dC;
    evec = vC;
else
    [evec, eval] = eig(C);
    deig         = diag(eval);
end

imag_e = find(abs(imag(deig))>0);
neg_e  = find(deig<0);
rem    = unique([imag_e;neg_e]);

if (length(rem)>0), 
  %warning('Covariance Matrix is not OK, redefined to be positive definite');
  deig(rem)   = [];
  evec(:,rem) = [];
end
k = length(deig);

proj = randn(N, k)*diag(sqrt(deig));
x    = ones(N, 1)*m + proj*evec';
