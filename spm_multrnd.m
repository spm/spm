function [m] = spm_multrnd(p,N)
% Sample from multinomial distribution
% FORMAT [m] = spm_multrnd(p,N)
%
% p    - [M x 1] vector of probabilities
% N    - Number of samples to generate
% 
% m    - [N x 1] vector of samples, where each sample is number from 1 to M
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_multrnd.m 3190 2009-06-08 17:13:36Z guillaume $

cp = [0; cumsum(p(:))];
m  = zeros(N,1);
for n=1:N
    m(n) = find(rand > cp, 1, 'last');
end
