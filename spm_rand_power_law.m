function [y] = spm_rand_power_law(G,Hz,dt,N)
% generates random variates with a power law spectral density
% FORMAT [y] = spm_rand_power_law(G,Hz,dt,N)
% G   - spectral densities (one per row)
% Hz  - frequencies
% dt  - sampling interval
% N   - number of time bins
%
% see also: spm_rand_mar
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_rand_power_law.m 5837 2014-01-18 18:38:07Z karl $
 

% create random process
%--------------------------------------------------------------------------
[m n] = size(G);
y     = randn(N,n)*sqrt(m);
s     = fft(y);
w     = (0:(N - 1))/dt/N;
g     = zeros(N,n);
for i = 1:m
    j      = find(w > Hz(i),1);
    g(j,:) = G(i,:);
end
y     = real(ifft(s.*sqrt(g)));



 
