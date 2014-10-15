function RCM = spm_sparse_regression(y,X,X0)
% Sparse (logistic) regression using Bayesian model reduction
% FORMAT RCM = spm_sparse_regression(y,X,X0)
% y   - univariate response variable
% X   - design matrix of explanatory variables
% X0  - confounds
% 
% RCM - reduced causal model structure
%   RCM.M      - GLM
%   RCM.Pp     - Model posterior (with and without each parameter)
%   RCM.Ep     - Bayesian parameter mean under reduced model
%   RCM.Cp     - Bayesian parameter covariance under reduced model
%   RCM.Vp     - Bayesian parameter variance under selected model
%__________________________________________________________________________
%
% spm_sparse_regression performs a sparse regression using priors on the
% parameters of a GLM and hyperpriors on the noise precision to recover a
% sparse set of explanatory variables. The implicit Bayesian model
% reduction (i.e., elimination of redundant parameters) uses post-hoc
% optimisation. If the response variable is in the range [0 1] then a logit
% transform is applied to produce sparse logistic regression.
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_sparse_regression.m 6247 2014-10-15 13:58:56Z guillaume $


% Logit transform
%--------------------------------------------------------------------------
if max(y) <= 1 && min(y) >= 0
    TOL = exp(-16);
    y   = log((y + TOL)./(1 - y + TOL));
end

% Scale
%--------------------------------------------------------------------------
scale  = std(y);
y      = y/scale;
X      = X/scale;

% Model specification
%--------------------------------------------------------------------------
M.IS   = @(P,M,U) U*P.A;
M.pE.A = zeros(size(X,2),1);
M.pC.A = M.pE.A + 1/mean(var(X))/16;
M.hE   = 4;
M.hC   = 1/64;
 
% Confounds (at the between subject level)
%--------------------------------------------------------------------------
try, Y.X0 = X0; end
Y.y       = y;
    
% Model inversion for this region
%--------------------------------------------------------------------------
[Ep,Cp] = spm_nlsi_GN(M,X,Y);
 
% Save model for subsequent BMS
%--------------------------------------------------------------------------
DCM.M   = M;
DCM.Ep  = Ep;
DCM.Cp  = Cp;
 

% Bayesian model reduction
%==========================================================================
RCM     = spm_dcm_post_hoc(DCM);

% Reorganise output structures
%--------------------------------------------------------------------------
RCM.Ep  = RCM.Ep.A;
RCM.Pp  = RCM.Pp.A;
RCM.Vp  = RCM.Vp.A;
RCM     = rmfield(RCM,'qP');
