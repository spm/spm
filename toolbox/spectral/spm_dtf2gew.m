function [gew,pve] = spm_dtf2gew(dtf,C)
% Converts directed transfer function to Geweke Granger causality
% FORMAT [gew,pve] = spm_csd2gew(dtf,C)
%
% dtf  (N,n,n)   - (unnormalised) directed or modulation functions
% C              - optional noise (fluctation) covariance matrix C(n,n)
%                - or cross spectral density C(N,n,n)
%                - or spectral power C(N,n)
%
% gwe  (N,n,n)   - Geweke's frequency domain Granger causality
% pve  (N,n,n)   - proportion of variance explained
%
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dtf2gew.m 6481 2015-06-16 17:01:47Z karl $


% preliminaries
%--------------------------------------------------------------------------
n  = size(dtf,2);
nw = size(dtf,1);

%  spectral density of fluctuations
%--------------------------------------------------------------------------
c     = zeros(nw,n,n);
if nargin < 2;
    C = eye(n,n);
end
if size(C,1) == nw
    for i = 1:n
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
    mtf        = squeeze(dtf(i,:,:));
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
        hkj        = abs(dtf(:,k,j)).^2;
        pve(:,k,j) = rkj.*hkj./sk;
        gew(:,k,j) = -log(1 - pve(:,k,j));
    end
end

