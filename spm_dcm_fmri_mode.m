function [U,E,F] = spm_dcm_fmri_mode(Ev,modes)
% Generates modes and matrices for spectral DCM from Lyapunov exponents
% FORMAT [Ep,Cp] = spm_dcm_fmri_mode_gen(Ev,modes,Cv)
% Ev    - (log of negative) Lyapunov exponents or eigenvalues of Jacobian
% modes - modes or eigenvectors
%
% U     - weighted modes; such that U*U' = F
% E     - (neuronal) effective  connectivity matrix
% F     - (neuronal) functional connectivity matrix E = -inv(F)/2
%
% This routine computes the connecivity graph for spectral DCM (modes).
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_mode.m 5823 2014-01-02 14:01:10Z guillaume $

    
% outer product
%==========================================================================
U    = modes*diag(sqrt(exp(Ev)/2));
E    = modes*diag(-exp(-Ev))*modes';
F    = modes*diag(exp(Ev)/2)*modes';
