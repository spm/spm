function DCM = spm_dcm_erp_update(DCM,oldDCM,fields)
% sets priors over DCM model parameters for Bayesian updating
% FORMAT DCM = spm_dcm_erp_update(DCM,oldDCM,fields)
%
% DCM    - DCM structure to be updated
% oldDCM - inverted DCM with posterior moments
% fields - character array of fields to be updated: e.g.,{'A','B'}
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_erp_sensitivity.m 4814 2012-07-30 19:56:05Z karl $

% get prior structure for the sort of model
%==========================================================================
try, model = DCM.options.model; catch, model = 'CMC'; end

% get the posterior (expectation and covariance) moments of parameters
%--------------------------------------------------------------------------
Ep      = oldDCM.Ep;
Cp      = spm_unvec(diag(oldDCM.Cp),Ep);
 
% prior moments of parameters
%--------------------------------------------------------------------------
[pE,pC] = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,model);

% fill-in prior moments
%--------------------------------------------------------------------------
for i = 1:length(fields)
    
    pE = setfield(pE,fields{i},getfield(Ep,fields{i}));
    pC = setfield(pC,fields{i},getfield(Cp,fields{i}));
    
end

% they surprising model structure
%==========================================================================
DCM.M.pE = pE;
DCM.M.pC = pC;