function [y] = spm_phase_shuffle(x)
% phase-shuffling of a vector
% FORMAT [y] = spm_phase_shuffle(x)
%__________________________________________________________________________


% randomise phase
%--------------------------------------------------------------------------
n            = length(x);
s            = fft(x);
i            = 2:ceil(n/2);
r            = rand(length(i),1)*2*pi - pi;
p            = zeros(n,1);
p(i)         =  r;
p(n - i + 2) = -r;
s            = abs(s).*exp(j*p);
y            = real(ifft(s));
