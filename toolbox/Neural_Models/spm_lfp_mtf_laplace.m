function [G,f] = spm_lfp_mtf_laplace(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum) 
% FORMAT [G,f] = spm_lfp_mtf(P,M,U)
%
% P - parameters
% M - neural mass model structure
% f - frequencies (per sec.)
%
% G - G(N,nc,nc} - cross-spectral density for nc channels
%                - for N frequencies in M.Hz [default 1:64Hz]
%                  
%__________________________________________________________________________
% Copyright (C) 2007 Wellcme Department of Imaging Neuroscienceo

% Rosalyn Moran
% $Id$
 
% compute log-spectral density
%==========================================================================
 
% frequencies of interest
%--------------------------------------------------------------------------
try
    dt = 1/(2*M.Hz(end));
    N  = 2*length(M.Hz(:));
catch
    N  = 128;
    dt = 1/N;
end
f    = [1:N/2]'/(N*dt);
w    = f*2*pi;
 
% spectrum of innovations or noise (U)
%--------------------------------------------------------------------------
U    = exp(P.a)*f.^(-1)*2;                 % spectral density of (AR) noise      
U    = U + exp(P.b);                       % spectral density of IID noise
    
% augment and bi-linearise
%--------------------------------------------------------------------------
[M0,M1,L] = spm_bireduce(M,P);
 
% compute modulation transfer function using FFT of the kernels
%--------------------------------------------------------------------------
[K0,K1]   = spm_kernels(M0,M1,L,N,dt);
 
% [cross]-spectral density
%--------------------------------------------------------------------------
[N,nc,nu] = size(K1);
nc    = length(P.L);
G     = zeros(N/2,nc,nc);

    A = full(M0(2:end,2:end)); 
    for i = 1:nc
        B(:,i) = full(M1{i}(2:end,1));
    end

    C = full(L(:,2:end));
    D = zeros(nc,nu); 

    for i = 1:nc
         for j = 1:nc
              for k =1:nu
                 [bi,ai]     =  spm_ss2tf(A,B,C(i,:),D(i,:),k); 
                 [bj,aj]     =  spm_ss2tf(A,B,C(j,:),D(j,:),k); 
                 Si          =  spm_freqs(bi,ai,w);
                 Sj          =  spm_freqs(bj,aj,w);
                 Gij         =  Si.*conj(Sj);
                 G(:,i,j)    =  G(:,i,j) + abs(Gij).*U*P.L(i,i)*P.L(j,j);

              end
         end
    end
