function [evidence] = spm_dcm_evidence (DCM)
% Compute evidence of DCM model
% FORMAT [evidence] = spm_dcm_evidence (DCM)
% 
% DCM       DCM data structure
%
% evidence   Contains the following fields:
%
%            .region_cost(i)        The cost of prediction errors in region i
%            .bic_penalty           Bayesian information criterion penalty 
%            .bic_overall           The overall BIC value
%            .aic_penalty           Akaike's information criterion penalty 
%            .aic_overall           The overall AIC value
%
%            All of the above are in units of NATS (not bits)
%
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id: spm_dcm_evidence.m 112 2005-05-04 18:20:52Z john $


v=DCM.v;

% Only look at those parameters with non-zero prior covariance
pCdiag=diag(DCM.M.pC);
wsel=find(~(pCdiag==0));

% Ensure residuals have been estimated correctly
if length(isnan(DCM.R))>0,
    % Assume that reason for incorrect R was zero columns in X0
    
    % 1. Remove zero columns from X0 if there are any 
    X0=DCM.Y.X0;
    ncol_X0=size(X0,2);
    new_X0=[];
    for col_X0=1:ncol_X0,
        if ~(length(find(X0(:,col_X0)==0))==DCM.v)
            new_X0=[new_X0 X0(:,col_X0)];
        end
    end
    X0=new_X0;
    
    % 2. Recompute residuals
    R     = DCM.Y.y - DCM.y;
    R     = R - X0*inv(X0'*X0)*(X0'*R);
    DCM.R = R;
    
    % Note DCM.R not saved to file
end

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





