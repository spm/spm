function [mtf,Hz] = spm_ssm2mtf(dfdx,dfdu,dgdx,Hz)
% computes cross spectral density from state space representation
% FORMAT [mtf,Hz] = spm_ssm2mtf(dfdx,dfdu,dgdx,Hz)
%
% dfdx - Jacobian
% dfdu - input matrix  [default: 1]
% dgdx - output matrix [default: 1]
% Hz   - frequencies   [default: based on maximum eigenvalue]
%
% mtf  - directed or modulation transfer function
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
if nargin < 2, dfdu = [];      end
if nargin < 3, dgdx = [];      end
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
s = 1j*imag(s) + min(real(s),-1/16);

% frequencies if unspecified
%--------------------------------------------------------------------------
if nargin < 4
    Hz  = -16*max(real(s));
    Hz  = Hz*(1:128)'/128;
end

% Modulaton transfer functions
%==========================================================================
nw    = size(Hz(:),1);        % number of frequencies
ng    = size(dgdx,1);         % number of outputs
nu    = size(dfdu,2);         % number of inputs
nk    = size(v,2);            % number of modes
mtf   = zeros(nw,ng,nu);

% derivatives over modes
%--------------------------------------------------------------------------
dgdv  = dgdx*v;
dvdu  = pinv(v)*dfdu;
for j = 1:nu
    for i = 1:ng
        for k = 1:nk
            
            % transfer functions (FFT of kernel)
            %--------------------------------------------------------------
            Sk         = 1./(1j*2*pi*Hz(:) - s(k));
            mtf(:,i,j) = mtf(:,i,j) + dgdv(i,k)*dvdu(k,j)*Sk;

        end
    end
end

