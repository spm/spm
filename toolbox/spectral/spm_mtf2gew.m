function [gew,pve] = spm_mtf2gew(mtf,C)
% Converts directed transfer function to Geweke Granger causality
% FORMAT [gew,pve] = spm_csd2gew(mtf,C)
%
% mtf  (N,n,n)   - (unnormalised) directed or modulation transfer function
% C              - optional noise (fluctation) covariance matrix C(n,n)
%                - or cross spectral density C(N,n,n)
%                - or spectral power C(N,n) [default: C = eye(n,n)]
%
% gew  (N,n,n)   - Geweke's frequency domain Granger causality
% pve  (N,n,n)   - proportion of variance explained
%
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
ns = size(mtf,3);
n  = size(mtf,2);
nw = size(mtf,1);

%  spectral density of fluctuations
%--------------------------------------------------------------------------
c     = zeros(nw,ns,ns);
if nargin < 2
    C = eye(ns,ns);
end
if size(C,1) == nw
    for i = 1:ns
        c(:,i,i) = C(:,i);
    end
end
if size(C,1) < nw
    for i = 1:nw
        c(i,:,:) = C;
    end
end

% cross spectral density
%--------------------------------------------------------------------------
csd   = zeros(nw,n,n);
for i = 1:nw
    mtf        = squeeze(mtf(i,:,:));
    C          = squeeze(c(i,:,:));
    csd(i,:,:) = mtf*C*mtf';
end

% Geweke Granger Causality in the Frequency domain
%--------------------------------------------------------------------------
pve   = zeros(nw,n,n);
gew   = zeros(nw,n,n);
for j = 1:n
    for k = 1:n
        rkj        = abs(c(:,j,j) - (c(:,j,k).^2)./c(:,k,k));
        sk         = abs(csd(:,k,k));
        hkj        = abs(mtf(:,k,j)).^2;
        pve(:,k,j) = rkj.*hkj./sk;
        gew(:,k,j) = -log(1 - pve(:,k,j));
    end
end
