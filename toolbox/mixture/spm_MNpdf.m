function [y] = spm_MNpdf (m, C, x)
% Evaluate a Multivariate Gaussian PDF
% FORMAT [y] = spm_MNpdf (m, C, x)
% 
% m     [d x 1] mean
% C     [d x d] covar
% x     [n x d] points at which to evaluate
%
% y     [n x 1] density at n points
%__________________________________________________________________________

% Will Penny 
% Copyright (C) 2007-2022 Wellcome Centre for Human Neuroimaging


ic     = inv(C);
[n, d] = size(x);

x    = x - ones(n, 1)*m(:)';
fact = sum(((x*ic).*x), 2);

y    = exp(-0.5*fact);
y    = y./sqrt((2*pi)^d*det(C));
