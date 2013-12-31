function [U,E,F] = spm_dcm_fmri_mode(Ev,modes)
% Generates modes and matrices for spectral DCM from Lyapunov exponents
% FORMAT [Ep,Cp] = spm_dcm_fmri_mode_gen(Ev,modes,Cv)
% Ev    - (log of negative) Lyapunov exponents or eigenvalues of Jacobian
% modes - modes or eigenvectors
%
% U     - wweighted modes; such that U*U' = F
% E     - (neuronal) effective  connectivity matrix
% F     - (neuronal) functional connectivity matrix E = -inv(F)/2
%
% This routine computes the connecivity graph for spectral DCM (modes)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_graph_gen.m 5817 2013-12-23 19:01:36Z karl $

    
% outer product
%==========================================================================
U    = modes*diag(sqrt(exp(Ev)/2));
E    = modes*diag(-exp(-Ev))*modes';
F    = modes*diag(exp(Ev)/2)*modes';


