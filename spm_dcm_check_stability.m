function [is_stable,eigval] = spm_dcm_check_stability(DCM)
% Checks the stability of a DCM model by examining the eigenvalue 
% spectrum for the intrinsic connectivity matrix (Lyapunov exponent)
%
% FORMAT spm_dcm_check_stability(DCM)
%
% DCM          model to check
% is_stable    returns 1 if stable, 0 if not stable
% eigval       lyapunov exponent
% 
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging
%
% Will Penny & Klaas Enno Stephan & Peter Zeidman

% Translate A-matrix to that used in spm_fx_fmri()
SE     = diag(DCM.Ep.A);
EE     = DCM.Ep.A - diag(exp(SE)/2 + SE);

% Check stability
eigval = max(eig(full(EE)));
is_stable = eigval < 0;

end