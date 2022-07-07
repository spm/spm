function [Ep,Cp] = spm_dcm_fmri_mode_gen(Ev,modes,Cv)
% Generate adjacency matrix for spectral DCM from Lyapunov exponents
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

% Karl Friston
% Copyright (C) 2013-2022 Wellcome Centre for Human Neuroimaging

    
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
