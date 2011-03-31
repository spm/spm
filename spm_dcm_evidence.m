function evidence = spm_dcm_evidence(DCM)
% Compute evidence of DCM model
% FORMAT evidence = spm_dcm_evidence(DCM)
%
% DCM       - DCM data structure
%
% evidence  - structure with the following fields
%   .region_cost(i)  - The cost of prediction errors in region i
%   .bic_penalty     - Bayesian information criterion penalty
%   .bic_overall     - The overall BIC value
%   .aic_penalty     - Akaike's information criterion penalty
%   .aic_overall     - The overall AIC value
%
% All of the above are in units of NATS (not bits).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_dcm_evidence.m 4281 2011-03-31 19:49:57Z karl $


% Only look at those parameters with non-zero prior covariance
%--------------------------------------------------------------------------
v     = DCM.v;                                   % number of samples
n     = DCM.n;                                   % number of regions
wsel  = find(diag(DCM.M.pC));

% Look at costs of coding prediction errors by region
%--------------------------------------------------------------------------
for i = 1:n
    try
        lambda_i = DCM.Ce(i*v,i*v);              % Ce is error covariance
    catch
        try
            lambda_i = DCM.Ce(i);                % Ce is a hyperparameter
        catch
            lambda_i = DCM.Ce;                   % Ce is the hyperparameter
        end
    end
    evidence.region_cost(i) = -0.5*v*log(lambda_i) ...
        - 0.5*DCM.R(:,i)'*(1/lambda_i)*eye(v)*DCM.R(:,i);
end

% Results
%--------------------------------------------------------------------------
evidence.aic_penalty = length(wsel);
evidence.bic_penalty = 0.5*length(wsel)*log(v);
evidence.aic_overall = sum(evidence.region_cost) - evidence.aic_penalty;
evidence.bic_overall = sum(evidence.region_cost) - evidence.bic_penalty;
