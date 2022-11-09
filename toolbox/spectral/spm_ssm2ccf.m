function [ccf,pst] = spm_ssm2ccf(dfdx,dfdu,dgdx,Hz)
% computes cross covariance from state space representation
% FORMAT [ccf,pst] = spm_ssm2ccf(dfdx,dfdu,dgdx,Hz)
%
% dfdx - Jacobian
% dfdu - input matrix  [default: 1]
% dgdx - output matrix [default: 1]
% Hz   - frequencies
%
% ccf  - cross covariance functions
% pst  - vector of lags for evaluation (seconds)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
if nargin < 2, dfdu = [];      end
if nargin < 3, dgdx = [];      end
if isempty(dfdu), dfdu = speye(size(dfdx)); end
if isempty(dgdx), dgdx = speye(size(dfdx)); end

% frequencies if unspecified
%--------------------------------------------------------------------------
if nargin < 4
    
    % condition unstable eigenmodes
    %----------------------------------------------------------------------
    s  = eig(full(dfdx),'nobalance');
    s  = min(real(s),-1/64);
    Hz = -16*max(real(s));
    Hz = Hz*(1:128)'/128;
    
end



% Cross covariance function by cross spectral density
%==========================================================================
[csd,Hz]  = spm_ssm2csd(dfdx,dfdu,dgdx,Hz);
[ccf,pst] = spm_csd2ccf(csd,Hz);