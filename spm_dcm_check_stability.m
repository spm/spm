function [is_stable,eigval] = spm_dcm_check_stability(DCM)
% Check stability of a DCM using Lyapunov exponent
% FORMAT [is_stable,eigval] = spm_dcm_check_stability(DCM)
%
% DCM        - DCM structure or its filename
%
% is_stable  - returns 1 if stable, 0 if not stable
% eigval     - Lyapunov exponent
%
% This function checks the stability of a DCM by examining the eigenvalue
% spectrum for the the intrinsic connectivity matrix (Lyapunov exponent).
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan & Peter Zeidman
% $Id: spm_dcm_check_stability.m 6077 2014-06-30 16:55:03Z spm $


%-Load DCM if necessary
%--------------------------------------------------------------------------
if ~isstruct(DCM)
    load(DCM);
end

%-Translate A-matrix to that used in spm_fx_fmri
%--------------------------------------------------------------------------
SE        = diag(DCM.Ep.A);
EE        = DCM.Ep.A - diag(exp(SE)/2 + SE);

%-Check stability
%--------------------------------------------------------------------------
eigval    = max(eig(full(EE)));
is_stable = eigval < 0;
