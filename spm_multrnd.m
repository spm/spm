function [m] = spm_multrnd(p,N)
% Sample from multinomial distribution
% FORMAT [m] = spm_multrnd(p,N)
%
% p    - [M x 1] vector of probabilities
% N    - Number of samples to generate
% 
% m    - [N x 1] vector of samples, where each sample is number from 1 to M
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2009-2022 Wellcome Centre for Human Neuroimaging


cp = [0; cumsum(p(:))];
m  = zeros(N,1);
for n=1:N
    m(n) = find(rand > cp, 1, 'last');
end
