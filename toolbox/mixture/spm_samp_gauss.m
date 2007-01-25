function [x] = spm_samp_gauss (m, C, N)
% Sample from a Gaussian PDF
% FORMAT [x] = spm_samp_gauss (m, C, N)
% m     [d x 1] mean
% C     [d x d] covar
% N     Number of samples
%
% x     [N x d] matrix of samples
%___________________________________________________________________________
% Copyright (C) 2007 Wellcome Department of Imaging Neuroscience

% Will Penny 
% $Id$

d = size(C, 1);
m = reshape(m, 1, d);   % Ensure that m is a row vector

[evec, eval] = eig(C);
deig=diag(eval);
if (~isreal(deig)) | any(deig<0), 
  warning('Covariance Matrix is not OK, redefined to be positive definite');
  eval=abs(eval);
end

proj = randn(N, d)*sqrt(eval);
x = ones(N, 1)*m + proj*evec';