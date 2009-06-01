function [m] = spm_multrnd(p,N)
% Sample from multinomial distribution
% FORMAT [m] = spm_multrnd(p,N)
%
% p     [M x 1] vector of probabilities
% N     Number of samples to generate
% 
% m     [N x 1] vector of samples, where each sample is number from 1 to M
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_multrnd.m 3170 2009-06-01 12:03:31Z will $

p=p(:);
p(end)=0;
p=circshift(p,1);
cp=cumsum(p);
for n=1:N,
    s=rand(1)>cp;
    s1=find(s==1);
    m(n)=max(s1);
end