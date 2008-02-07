function [x,label] = spm_samp_mix (mix, N)
% Sample from a Gaussian Mixture PDF
% FORMAT [x,label] = spm_samp_mix (mix, N)
%
% mix   Data structure for mixture model (see spm_mix for info)
% N     Number of samples
%
% x     [N x d] matrix of samples
% label [N x 1] vector of sample labels
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_samp_mix.m 1143 2008-02-07 19:33:33Z spm $

priors = rand(1, N);
cum_prior = 0;      
total_samples = 0;  
label = zeros(N, 1);
for j=1:mix.m,
    num_samples = sum(priors >= cum_prior & ...
        priors < cum_prior + mix.state(j).prior);
    m=mix.state(j).m;
    C=squeeze(mix.state(j).C);
    x(total_samples+1:total_samples+num_samples, :) = ...
        spm_samp_gauss(m, C, num_samples);
    cum_prior = cum_prior + mix.state(j).prior;
    total_samples = total_samples + num_samples;
    label(total_samples+1:total_samples+num_samples) = j;
end