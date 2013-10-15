function [y] = spm_dcm_fmri_csd_gen(x,v,P)
% Conversion routine for DEM inversion of DCM for CSD (fMRI)
% FORMAT [y] = spm_dcm_fmri_csd_gen(x,v,P)
%
% This routine computes the spectral response of a network of regions
% driven by  endogenous fluctuations and exogenous (experimental) inputs.
% It returns the complex cross spectra of regional responses as a
% three-dimensional array. The endogenous innovations or fluctuations are
% parameterised in terms of a (scale free) power law, in frequency space.
%
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_csd_gen.m 5696 2013-10-15 19:10:26Z karl $


% global DCM and evaluate orginal generative model
%==========================================================================
global DCM

% conditional prediction
%--------------------------------------------------------------------------
y   = feval(DCM.M.IS,v,DCM.M,DCM.U);
y   = feval(DCM.M.FS,y,DCM.M);
y   = spm_vec(y);
