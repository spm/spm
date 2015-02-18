function [post,exp_r,xp,pxp,bor] = spm_dcm_bmc(DCM)
% Bayesian model comparison
% FORMAT [post,exp_r,xp,pxp,bor] = spm_dcm_bmc(DCM)
%
% DCM     - {subjects x models} cell array of DCMs
% ------------------------------------------------
%     DCM{i,j}.F  - free energy
%
% OUTPUTS
% -------
% post    - FFX posterior model probabilities p(m|y)
% exp_r   - RFX expectation of the posterior  p(m|y)
% xp      - RFX exceedance probabilities
% pxp     - RFX protected exceedance probabilities
% bor     - RFX Bayes Omnibus Risk (probability that model frequencies 
%           are equal)
%
% This routine computes fixed and random effects posterior probabilities
% over models. It also returns exceedance  probabilities and protected
% statistics.
% 
% See also: spm_dcm_bma.m and spm_BMS.m
%__________________________________________________________________________
% 
% References:
%
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009)
% Bayesian Model Selection for Group Studies. NeuroImage 46:1004-1017
%
% Rigoux, L, Stephan, KE, Friston, KJ and Daunizeau, J. (2014)
% Bayesian model selection for group studies - Revisited. 
% NeuroImage 84:971-85. doi: 10.1016/j.neuroimage.2013.08.065
%__________________________________________________________________________
% Copyright (C) 2009-2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_dcm_bmc.m 6343 2015-02-18 16:46:00Z spm $


% assemble log evidence
%--------------------------------------------------------------------------
[n,m] = size(DCM);
for i = 1:n
    for j = 1:m
        F(i,j) = DCM{i,j}.F;
    end
end

% FFX posterior over models
%--------------------------------------------------------------------------
P    = sum(F,1);
P    = P - max(P);
P    = exp(P);
post = P/sum(P);

% RFX posterior over models
%--------------------------------------------------------------------------
if nargout > 1
    [alpha,exp_r,xp,pxp,bor] = spm_BMS(F);
end
