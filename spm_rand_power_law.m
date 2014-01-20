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
% $Id: spm_rand_power_law.m 5841 2014-01-20 10:19:25Z karl $
 

% create random process
%--------------------------------------------------------------------------
[m n] = size(G);
w     = (0:(N - 1))/dt/N;
dHz   = Hz(2) - Hz(1);
g     = zeros(N,n);
for i = 1:m
    j      = find(w > Hz(i),1);
    s      = sqrt(G(i,:)).*(randn(1,n) + 1j*randn(1,n));
    g(j,:) = s;
    j      = N - j + 2;
    g(j,:) = conj(s);
    
end
y     = real(ifft(g));
y     = y*sqrt(mean(sum(G)*dHz)/mean(var(y)));
