function [evidence] = spm_dcm_evidence (DCM)

% function [evidence] = spm_dcm_evidence (DCM)
%
% Compute evidence of DCM model
% 
% DCM       DCM data structure
%
% evidence   Contains the following fields:
%
%            .region_cost(i)        The cost of prediction errors in region i
%            .bic_penalty           Bayesian information criterion penalty (in units of NATS)
%            .bic_overall           The overall BIC value
%            .aic_penalty           Akaike's information criterion penalty (in units of NATS)
%            .aic_overall           The overall AIC value
%
%            All of the above are in units of NATS (not bits)
%
% %W% Will Penny %E%

v=DCM.v;

% Only look at those parameters with non-zero prior covariance
pCdiag=diag(DCM.M.pC);
wsel=find(~(pCdiag==0));


% Look at costs of coding prediction errors by region
n=size(DCM.A,1); 
for i=1:n,
    lambda_i=DCM.Ce(i*v,i*v);
    evidence.region_cost(i)=-0.5*v*log(lambda_i);
    evidence.region_cost(i)=evidence.region_cost(i)-0.5*DCM.R(:,i)'*(1/lambda_i)*eye(v)*DCM.R(:,i);
end

evidence.aic_penalty=length(wsel);
evidence.bic_penalty=0.5*length(wsel)*log(v);

evidence.aic_overall=sum(evidence.region_cost)-evidence.aic_penalty;
evidence.bic_overall=sum(evidence.region_cost)-evidence.bic_penalty;





