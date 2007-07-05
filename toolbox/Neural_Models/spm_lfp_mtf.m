function [G,f] = spm_lfp_mtf(P,M,U)
% Spectral response of a NMM (Modulation transfer function plus noise) 
% FORMAT [G,f] = spm_lfp_mtf(P,M,U)
%
% P - parameters
% M - neural mass model stucture
% f - ferquencies (per sec.)
%
% G - G(f); f = 1:64 Hz
%__________________________________________________________________________


% compute log-spectral density
%==========================================================================
i     = 1;    % i-th input
j     = 1;    % j-th output

% frequencies of interest
%--------------------------------------------------------------------------
try
    f    = M.Hz;
    dt   = 1/M.fs;
catch
    N    = 128;
    dt   = 1/N;
    f    = [1:N/2]/(N*dt);
end

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
G         = abs(S).^2;                             % energy from neural mass
warning on

% spectral density of (AR) noise
%--------------------------------------------------------------------------
try                      
    G   = G + diag(P.L*P.L')*exp(P.a)*f.^(-1)/64;  % energy from 1/f process
end

% spectral density of IID noise
%--------------------------------------------------------------------------
try
    G   = G + diag(P.L*P.L')*exp(P.b)*(f.^0)/1024; % energy from IID process
end


return

% NB: compute modulation transfer function using FFT of the kernels
%--------------------------------------------------------------------------
[K0,K1] = spm_kernels(M0,M1,L1,N,dt);
S       = fft(K1(:,:,i));
G       = abs(S([1:N/2] + 1,j)).^2;



