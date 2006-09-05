function [G,f] = spm_lfp_mtf(P,M)
% Modulation transfer function g(f) of a NMM
% FORMAT [G,f] = spm_lfp_mtf(P,M)
%
% P - parameters
% M - neural mass model stucture
% f - ferquencies (per sec.)
%
% G - G(f); f = 1:64 Hz
%__________________________________________________________________________


% compute log-spectral density
%==========================================================================
i       = 1;    % i-th input
j       = 1;    % j-th output

% frequnecies of interest
%--------------------------------------------------------------------------
N       = 128;
dt      = 1/N;
f       = [1:N/2]/(N*dt);

% augment and bi-linearise
%--------------------------------------------------------------------------
[M0,M1,L1] = spm_bireduce(M,P);

% compute modulation transfer function using FFT of the kernels
%--------------------------------------------------------------------------
warning off
A         = full(M0(2:end,2:end));
B         = full(M1{i}(2:end,1));
C         = full(L1(j,2:end));
[b,a]     = ss2tf(A,B,C,0,i);
[num,den] = bilinear(b,a,1/dt);
[S f]     = freqz(num,den,f,1/dt);
G         = abs(S).^2;
warning on

return

% compute modulation transfer function using FFT of the kernels
%--------------------------------------------------------------------------
[K0,K1] = spm_kernels(M0,M1,L1,N,dt);
S       = fft(K1(:,:,i));
G       = abs(S([1:N/2] + 1,j)).^2;



