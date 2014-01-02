function [Ep,Cp] = spm_dcm_fmri_mode_gen(Ev,modes,Cv)
% Generates adjacency matrix for spectral DCM from Lyapunov exponents
% FORMAT [Ep,Cp] = spm_dcm_fmri_mode_gen(Ev,modes,Cv)
% Ev    - Lyapunov exponents or eigenvalues of effective connectivity
% modes - modes or eigenvectors
% Cv    - optional (posterior) covariance matrix
%
% Ep    - Jacobian or (symmetric) effective connectivity matrix
% Cp    - posterior covariance matrix of Jacobian elements
%
% This routine computes the connecivity graph for spectral DCM (modes).
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_mode_gen.m 5823 2014-01-02 14:01:10Z guillaume $

    
% outer product
%==========================================================================
Ep    = modes*diag(-exp(-Ev))*modes';

if nargout == 1, return, end

% covariance
%==========================================================================
dAdv  = spm_diff(@spm_dcm_fmri_mode_gen,Ev,modes,1);
for i = 1:length(dAdv)
    G(:,i) = dAdv{i}(:);
end
Cp    = G*Cv*G';
