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
% model
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_eeg_channelmodes.m 4768 2012-06-11 17:06:55Z karl $
 
% number of channels and modes
%--------------------------------------------------------------------------
if nargin < 2, Nm = 8; end
 
% Spatial modes
%--------------------------------------------------------------------------
Nc    = dipfit.Nc;
pE    = spm_L_priors(dipfit);
if Nc < Nm
    
    % number of channels is already les than the number of modes
    %----------------------------------------------------------------------
    U     = speye(Nc,Nc);
else
    
    % evaluate eigenmodes of gain of covariance in sensor space
    %----------------------------------------------------------------------
    dGdg  = spm_diff('spm_lx_erp',pE,dipfit,1);
    L     = spm_cat(dGdg);
    U     = spm_svd(L*L',exp(-8));
    try
        U = U(:,1:Nm);
    end
end
