function y = spm_dcm_fmri_csd_gen(x,v,P)
% Conversion routine for DEM inversion of DCM for CSD (fMRI)
% FORMAT y = spm_dcm_fmri_csd_gen(x,v,P)
%
% This routine computes the spectral response of a network of regions
% driven by  endogenous fluctuations and exogenous (experimental) inputs.
% It returns the complex cross spectra of regional responses as a
% three-dimensional array. The endogenous innovations or fluctuations are
% parameterised in terms of a (scale free) power law, in frequency space.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2013-2022 Wellcome Centre for Human Neuroimaging


% global DCM and evaluate original generative model
%==========================================================================
global GLOBAL_DCM

% conditional prediction
%--------------------------------------------------------------------------
y = feval(GLOBAL_DCM.M.IS,v,GLOBAL_DCM.M,GLOBAL_DCM.U);
y = feval(GLOBAL_DCM.M.FS,y,GLOBAL_DCM.M);
y = spm_vec(y);
