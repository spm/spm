function DCM = spm_dcm_erp_update(DCM,oldDCM,fields)
% Set priors over DCM model parameters for Bayesian updating
% FORMAT DCM = spm_dcm_erp_update(DCM,oldDCM,fields)
%
% DCM    - DCM structure to be updated
% oldDCM - inverted DCM with posterior moments
% fields - character array of fields to be updated: e.g.,{'A','B'}
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_erp_update.m 5801 2013-12-10 18:45:20Z guillaume $

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
    
    pE.(fields{i}) = Ep.(fields{i});
    pC.(fields{i}) = Cp.(fields{i});
    
end

% they surprising model structure
%==========================================================================
DCM.M.pE = pE;
DCM.M.pC = pC;
