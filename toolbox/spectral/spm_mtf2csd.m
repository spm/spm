function [csd] = spm_mtf2csd(mtf,C)
% Converts modulation transfer function to cross spectral density
% FORMAT [csd] = spm_mtf2csd(mtf,C)
%
% mtf  (N,n,n)   - (unnormalised) directed or modulation transfer function
% C              - optional noise (fluctation) covariance matrix C(n,n)
%                - or cross spectral density C(N,n,n)
%                - or spectral power C(N,n) [default: C = eye(n,n)]
%
% csd  (N,n,n)   - cross spectral density
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
else
    for i = 1:nw
        c(i,:,:) = C(1:ns,1:ns);
    end
end

% cross spectral density
%--------------------------------------------------------------------------
csd   = zeros(nw,n,n);
for i = 1:nw
    tf         = squeeze(mtf(i,:,:));
    C          = squeeze(c(i,:,:));
    csd(i,:,:) = tf*C*tf';
end
