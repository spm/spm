function [gew,pve,H] = spm_csd2gew(csd,Hz,u)
% Convert cross sspectral density to Geweke Granger causality
% FORMAT [gew,pve,H] = spm_csd2gew(csd,Hz)
%
% ccf  (N,m,m)   - cross covariance functions
% Hz   (n x 1)   - vector of frequencies (Hz)
% u    (1)       - regularizer (default: 1);
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
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd2gew.m 5908 2014-03-05 20:31:57Z karl $

% preliminaries
%--------------------------------------------------------------------------
if nargin < 3
    try
        [gew,pve,H] = spm_csd2gew(csd,Hz,2);
        return
    catch
        [gew,pve,H] = spm_csd2gew(csd,Hz,16);
        return
    end
end

% pad spectrum if necessary
%--------------------------------------------------------------------------
iw    = 1 + round(Hz/(Hz(2) - Hz(1)));
n     = size(csd,2);
nw    = 257;
is    = ceil(nw/2):nw;

% Wilson-Burg algorithm
%==========================================================================

% initialise transfer function
%--------------------------------------------------------------------------
H     = zeros(nw,n,n);
P     = zeros(nw,n,n);
e     = norm(squeeze(max(csd,[],1)))/128;
E     = eye(n,n)*e;

P(iw,:,:) = csd;
for i = 1:n
    P(:,i,i) = P(:,i,i) + e;
    H(:,i,i) = sqrt(P(:,i,i));
end

% iterate until convergence: solve for H*H' = P
%--------------------------------------------------------------------------
for t = 1:128
    
    % compute left-hand side (deconvolution)
    %----------------------------------------------------------------------
    for w = 1:nw
        S        = squeeze(H(w,:,:)) + E;
        A(w,:,:) = (eye(n,n) + S*S'\squeeze(P(w,:,:)));
    end
    
    % retain causal signal and half zero lag
    %----------------------------------------------------------------------
    S           = ifft(A);
    S(is,:,:)   = 0;
    S(1,:,:)    = S(1,:,:)/2;
    A           = fft(S);
    
    % recover next update (convolution)
    %----------------------------------------------------------------------
    nrm   = zeros(nw,1);
    for w = 1:nw
        U        = squeeze(A(w,:,:));
        H(w,:,:) = squeeze(H(w,:,:))*U^(1/u);
        nrm(w)   = nrm(w) + norm(eye(n,n) - U,'inf');
    end
    
    % break if convergence
    %----------------------------------------------------------------------
    nrm = mean(nrm);
    if nrm < 1e-6, break, end
    if nrm > 8,   return, end
   
end

% transfer function and noise covariance
%==========================================================================

% get noise covariance
%--------------------------------------------------------------------------
S     = ifft(H);
R     = squeeze(S(1,:,:));
C     = real(R*R');
c     = sqrtm(C);

% recover transfer function
%--------------------------------------------------------------------------
for w = 1:nw
    H(w,:,:) = squeeze(H(w,:,:))/c;
    P(w,:,:) = squeeze(H(w,:,:))*C*squeeze(H(w,:,:))';
end

% Geweke Granger Causality in the Frequency domain
%--------------------------------------------------------------------------
pve   = zeros(nw,n,n);
gew   = zeros(nw,n,n);
for j = 1:n
    for k = 1:n
        rkj        = C(j,j) - (C(j,k)^2)/C(k,k);
        sk         = abs(P(:,k,k));
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
