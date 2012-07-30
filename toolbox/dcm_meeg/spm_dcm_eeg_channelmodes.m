function [U] = spm_dcm_eeg_channelmodes(dipfit,Nm)
% Returns the channel eigenmodes
% FORMAT [U] = spm_dcm_eeg_channelmodes(dipfit,Nm)
% dipfit  - spatial model specification
% Nm      - number of modes required (upper bound)
% U       - channel eigenmodes
%__________________________________________________________________________
%
% Uses SVD to identify the patterns with the greatest prior covariance
% assuming independent source activity in the specified spatial (forward)
% model. 
%
% U is scaled to ensure trace(U'*L*L'*U) = Nm
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_eeg_channelmodes.m 4814 2012-07-30 19:56:05Z karl $
 
% number of channels and modes
%--------------------------------------------------------------------------
if nargin < 2, Nm = 8; end

% Spatial modes
%--------------------------------------------------------------------------
pE    = spm_L_priors(dipfit);

% evaluate eigenmodes of gain of covariance in sensor space
%--------------------------------------------------------------------------
dGdg  = spm_diff('spm_erp_L',pE,dipfit,1);
L     = spm_cat(dGdg);
[U S] = spm_svd(L*L',exp(-8));
S     = diag(S);

% eigen-mode reduction
%--------------------------------------------------------------------------
try
    U = U(:,1:Nm);
    S = S(  1:Nm);
end

% re-scale spatial projector
%--------------------------------------------------------------------------
U     = U/sqrt(mean(S));

