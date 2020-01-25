function [ker,pst] = spm_ssm2ker(dfdx,dfdu,dgdx,pst)
% computes cross spectral density from state space representation
% FORMAT [ker,pst] = spm_ssm2ker(dfdx,dfdu,dgdx,pst)
%
% dfdx - Jacobian
% dfdu - input matrix  [default: 1]
% dgdx - output matrix [default: 1]
% pst  - time          [default: based on maximum eigenvalue]
%
% ker  - first-order (Volterra) kernels
%
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ssm2ker.m 7774 2020-01-25 18:07:03Z karl $

% preliminaries
%--------------------------------------------------------------------------
if nargin < 2,    dfdu = [];                end
if nargin < 3,    dgdx = [];                end
if isempty(dfdu), dfdu = speye(size(dfdx)); end
if isempty(dgdx), dgdx = speye(size(dfdx)); end

% Jacobian and eigenspectrum
%==========================================================================
try
   [v,s] = eig(full(dfdx),'nobalance');
catch
   v  = eye(size(dfdx));
   s  = NaN(size(dfdx));
end
s     = diag(s);

% condition unstable eigenmodes
%--------------------------------------------------------------------------
s = 1j*imag(s) + min(real(s),-1/64);

% time lags if unspecified
%--------------------------------------------------------------------------
if nargin < 4
    pst   = -8/max(real(s));
    pst   = pst*(1:64)'/64;
end

% Transfer functions
%==========================================================================

% transfer functions (FFT of kernel)
%--------------------------------------------------------------------------
ng    = size(dgdx,1);         % number of outputs
nu    = size(dfdu,2);         % number of inputs
nk    = size(v,2);            % number of modes
nt    = size(pst(:),1);       % number of time points
ker   = zeros(nt,ng,nk);      % kernels

% derivatives over modes
%--------------------------------------------------------------------------
dgdv  = dgdx*v;
dvdu  = pinv(v)*dfdu;
for j = 1:nu
    for i = 1:ng
        for k = 1:nk
            
            % kernels
            %--------------------------------------------------------------
            Kk         = exp(s(k)*pst(:));
            ker(:,i,j) = ker(:,i,j) + real(dgdv(i,k)*dvdu(k,j)*Kk);
            
        end
    end
end

