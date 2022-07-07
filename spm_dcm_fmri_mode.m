function [U,E,F] = spm_dcm_fmri_mode(Ev,modes)
% Generate modes and matrices for spectral DCM from Lyapunov exponents
% FORMAT [U,E,F] = spm_dcm_fmri_mode(Ev,modes)
% Ev    - (log of negative) Lyapunov exponents or eigenvalues of Jacobian
% modes - modes or eigenvectors
%
% U     - weighted modes; such that U*U' = F
% E     - (neuronal) effective  connectivity matrix
% F     - (neuronal) functional connectivity matrix E = -inv(F)/2
%
% This routine computes the connecivity graph for spectral DCM (modes).
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2013-2022 Wellcome Centre for Human Neuroimaging

    
% outer product
%==========================================================================
U    = modes*diag(sqrt(exp(Ev)/2));
E    = modes*diag(-exp(-Ev))*modes';
F    = modes*diag(exp(Ev)/2)*modes';
