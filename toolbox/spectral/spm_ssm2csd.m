function [csd,Hz] = spm_ssm2csd(dfdx,dfdu,dgdx,Hz)
% computes cross spectral density from state space representation
% FORMAT [csd,Hz] = spm_ssm2csd(dfdx,dfdu,dgdx,Hz)
%
% dfdx - Jacobian
% dfdu - input matrix  [default: 1]
% dgdx - output matrix [default: 1]
% Hz   - frequencies   [default: based on maximum eigenvalue]
%
% csd  - cross spectral density
%
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ssm2csd.m 7774 2020-01-25 18:07:03Z karl $

% preliminaries
%--------------------------------------------------------------------------
if nargin < 2,    dfdu = [];                end
if nargin < 3,    dgdx = [];                end
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


% cross spectral density via modulation transfer function
%==========================================================================
mtf = spm_ssm2mtf(dfdx,dfdu,dgdx,Hz);
csd = spm_mtf2csd(mtf);