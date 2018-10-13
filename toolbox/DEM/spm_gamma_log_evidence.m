function [F,sA] = spm_gamma_log_evidence(qA,pA,rA,~)
% Bayesian model reduction for gamma distibutions
% FORMAT [F,sA] = spm_gamma_log_evidence(qA,pA,rA,'shape')
%
% qA  - shape/rate parameter of posterior of full model
% pA  - shape/rate parameter of prior of full model
% rA  - shape/rate parameter of prior of reduced model
%
% 'shape' - for shape parameters (otherwise, rate parameters are assumed)
%
% F   - (negative) free energy or log evidence of reduced model
% sA  - shape/rate parameter of reduced posterior
%
% This routine compute the negative log evidence of a reduced model of a
% gamma distribution parameterised in terms of its shape parameter.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gamma_log_evidence.m 7447 2018-10-13 15:32:10Z karl $

if nargin < 4
    
    % change in free energy or log model evidence - rate parameter
    %--------------------------------------------------------------------------
    sA = qA + rA - pA;
    F  = - log(qA) - log(rA) + log(pA) + log(sA);
    
else
    
    % change in free energy or log model evidence - shape parameter
    %--------------------------------------------------------------------------
    sA = qA + rA - pA;
    F  = - gammaln(qA) - gammaln(rA) + gammaln(pA) + gammaln(sA);
    
end




