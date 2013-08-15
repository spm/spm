function [y] = spm_rand_power_law(m,n,a)
% Generate random variates with a power law spectral density
% FORMAT [y] = spm_rand_power_law(m,n,a)
% m   - time bins
% n   - varies
% a   - power law exponent g(w) = c*w^(-a): sum(g(w)) = 1
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_rand_power_law.m 5615 2013-08-15 14:37:24Z spm $
 

% create random process
%--------------------------------------------------------------------------
w              = (2:ceil(m/2))';
r              = rand(length(w),n)*2*pi - pi;
f              = w.^(-a/2);
f              = sqrt(m)*f/mean(sqrt(f));
q              = randn(length(w),n).*(f*ones(1,n));
p              = zeros(m,n);
p(w,:)         = q.*exp(1j*r);
p(m - w + 2,:) = q.*exp(-1j*r);
y              = real(ifft(p));
