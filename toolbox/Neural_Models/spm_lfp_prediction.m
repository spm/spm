function [G,w] = spm_lfp_prediction(P,M,varargin)
% prediction for log-spectral density of a NMM
% FORMAT [G,w] = spm_lfp_prediction(P,M)
%
% P - parameters
% M - neural mass model stucture
%
% G - ln(G(w))
%__________________________________________________________________________


% compute log-spectral density
%==========================================================================

% augment and bi-linearise
%--------------------------------------------------------------------------
[M0,M1,L1] = spm_bireduce(M,P);

% compute transfer functions
%--------------------------------------------------------------------------
N          = 128;
dt         = 8/1000;
t          = 1:(N/2);
w          = (t - 1)/(N*dt);
[K0,K1]    = spm_kernels(M0,M1,L1,N,dt);
S          = fft(K1);
G          = log((abs(S(t,:)).^2));



