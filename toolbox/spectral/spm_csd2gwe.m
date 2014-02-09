function [gew,pve,H] = spm_csd2gwe(csd,Hz)
% Converts cross sspectral density to Geweke Granger causality
% FORMAT [gew,pve,H] = spm_csd2gwe(csd,Hz)
%
% ccf  (N,m,m)   - cross covariance functions
% Hz   (n x 1)   - vector of frequencies (Hz)
%
% gwe  (N,m,m)   - Geweke's frequency domain Granger causality
% pve  (N,m,m)   - proportion of variance explained
% H    (N,m,m)   - transfer function matrix
%
% This routine uses the Wilson-Burg algorithm to perform spectral matrix
% factorisation. The minimum phase factor is then used to form the noise
% covariance (covariance of the innovations) and implicitly derive the
% transfer functions (and spectral Granger causality).
%
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ccf2gwe.m 5853 2014-01-24 20:38:11Z karl $

% pad spectrum if necessary
%--------------------------------------------------------------------------
iw          = 1 + round(Hz/(Hz(2) - Hz(1)));
csd(iw,:,:) = csd;
[nw,d,n]    = size(csd);
is          = round(nw/2):nw;

% Wilson-Burg algorithm
%==========================================================================

% initialise transfer function
%--------------------------------------------------------------------------
for w = 1:nw
    H(w,:,:) = eye(n,n);
end

% iiterate until convergence
%--------------------------------------------------------------------------
for t = 1:32
    
    % compute left-hand side (deconvolution)
    %----------------------------------------------------------------------
    for w = 1:nw
       S        = squeeze(H(w,:,:));
       A(w,:,:) = (eye(n,n) + S*S'\squeeze(csd(w,:,:)));
    end
    
    % retain causal signal and half zero lag
    %----------------------------------------------------------------------
    S           = ifft(A);
    S(is,:,:)   = 0;
    S(1,:,:)    = S(1,:,:)/2;
    A           = fft(S);
 
    % recover next update (convolution)
    %----------------------------------------------------------------------
    nrm   = 0;
    for w = 1:nw
       U        = squeeze(A(w,:,:));
       H(w,:,:) = squeeze(H(w,:,:))*U;
       nrm      = nrm + norm(eye(n,n) - U,'inf');
    end
    
    % break if convergence
    %----------------------------------------------------------------------
    if nrm < 1e-4, break, end
    
end

% transfer function and noise covariance
%==========================================================================

% get noise covariance
%--------------------------------------------------------------------------
S     = ifft(A);
R     = squeeze(S(1,:,:));
C     = real(R*R');

% recover transfer function
%--------------------------------------------------------------------------
for w = 1:nw
    H(w,:,:) = squeeze(H(w,:,:))/C;
end

% Geweke Granger Causality in the Frequency domain
%--------------------------------------------------------------------------
for j = 1:d
    for k = 1:d
        rkj        = C(j,j) - (C(j,k)^2)/C(k,k);
        sk         = abs(csd(:,k,k));
        hkj        = abs(H(:,k,j)).^2;
        pve(:,k,j) = rkj*hkj./sk;
        gew(:,k,j) = -log(1 - pve(:,k,j));
    end
end

% return  specified frequencies
%--------------------------------------------------------------------------
gew = gew(iw,:,:);
pve = pve(iw,:,:);
H   = H(iw,:,:);




