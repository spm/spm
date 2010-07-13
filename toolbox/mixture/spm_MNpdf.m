function [y] = spm_MNpdf (m, C, x)
% Evaluate a Multivariate Gaussian PDF
% FORMAT [y] = spm_MNpdf (m, C, x)
% 
% m     [d x 1] mean
% C     [d x d] covar
% x     [n x d] points at which to evaluate
%
% y     [n x 1] density at n points
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_MNpdf.m 3995 2010-07-13 17:19:49Z karl $

ic     = inv(C);
[n, d] = size(x);

x    = x - ones(n, 1)*m(:)';
fact = sum(((x*ic).*x), 2);

y    = exp(-0.5*fact);
y    = y./sqrt((2*pi)^d*det(C));
