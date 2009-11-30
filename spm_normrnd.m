function x = spm_normrnd(m, C, N)
% Random samples from Gaussian distribution 
% FORMAT x = spm_normrnd(m, C, N)
% m        - [d x 1] mean
% C        - [d x d] covariance or cell array {dC, vC} so that
%            [vC, diag(dC)] = eig(C)
% N        - number of samples
%
% x        - [d x N] matrix of samples
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_normrnd.m 3603 2009-11-30 18:56:50Z guillaume $

if iscell(C)
    deig = C{1};
    evec = C{2};
else
    [evec, eval] = eig(C);
    deig         = diag(eval);
end

i = (abs(imag(deig))>0) | (deig<0);
if any(i)
  %warning('Covariance matrix is not positive semi-definite: redefined');
  deig(i)   = [];
  evec(:,i) = [];
end

proj = randn(N,length(deig)) * diag(sqrt(deig));
x    = repmat(m(:),1,N) + evec*proj';
