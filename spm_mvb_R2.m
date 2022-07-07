function [R2,X,P] = spm_mvb_R2(MVB)
% Return the proportion of variance explained by the (MVB) MAP estimates
% FORMAT [R2,X,P] = spm_mvb_R2(MVB)
%
% MVB - MVB structure
% R2  - proportion of variance explained
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2011-2022 Wellcome Centre for Human Neuroimaging
 

% MAP predictions
%--------------------------------------------------------------------------
Y  = MVB.Y*MVB.M.qE;
 
% target variable
%--------------------------------------------------------------------------
X  = MVB.X;
 
% linearly optimised predictor
%--------------------------------------------------------------------------
Y  = [Y MVB.X0];
P  = Y*pinv(Y)*X;
 
% R-squared
%--------------------------------------------------------------------------
R2 = var(P)/var(X);
