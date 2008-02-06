function [y,f] = spm_lfp_mtf(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum) 
% FORMAT [G,f] = spm_lfp_mtf(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
%
% G - {G(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
%                  
%__________________________________________________________________________
 % Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_lfp_mtf.m 1132 2008-02-06 14:12:17Z karl $

 
% compute log-spectral density
%==========================================================================
 
% frequencies of interest
%--------------------------------------------------------------------------
try
    dt = 1/(2*M.Hz(end));
    N  = 1/dt;
    If = M.Hz(1):M.Hz(end);
catch
    N  = 128;
    dt = 1/N;
    If = 1:N/2;
end
f    = [1:N/2]';
 
 
% spectrum of innovations or noise (Gu)
%--------------------------------------------------------------------------
Gu   = exp(P.a)*f.^(-1)*2;                % spectral density of (AR) input      
Gu   = Gu + exp(P.b);                     % spectral density of IID input
 
% trial-specific effects
%==========================================================================
try, X = U.X; catch, X = sparse(1,0); end
for  c = 1:size(X,1)
    
    % exogenous and trial-specific inputs
    %----------------------------------------------------------------------
    M.u  = [sparse(length(P.C),1) X(c,:)];
 
    % augment and bi-linearise
    %----------------------------------------------------------------------
    [M0,M1,L] = spm_bireduce(M,P);
    
    % project onto spatial modes
    %----------------------------------------------------------------------
    try
        L = M.U'*L;
    end
 
    % compute modulation transfer function using FFT of the kernels
    %----------------------------------------------------------------------
    [K0,K1]   = spm_kernels(M0,M1,L,N,dt);
    
 
    % [cross]-spectral density
    %----------------------------------------------------------------------
    [N,nc,nu] = size(K1);
    G     = zeros(N/2,nc,nc);
    for i = 1:nc
        for j = i:nc
            
            % exogenous components
            %--------------------------------------------------------------
            for k = 1:nu
                Si       = fft(K1(:,i,k));
                Sj       = fft(K1(:,j,k));
                Gij      = Si.*conj(Sj);
                Gij      = abs(Gij([1:N/2] + 1)).*Gu;
                G(:,i,j) = G(:,i,j) + Gij;
                G(:,j,i) = G(:,j,i) + Gij;
            end
        end
    end
 
    % save frequencies of interest
    %----------------------------------------------------------------------
    y{c} = G(If,:,:);
    
end
 
% prediction
%--------------------------------------------------------------------------
if ~length(U)
    y = y{1};
end
