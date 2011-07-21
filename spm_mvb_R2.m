function [R2,X,P] = spm_mvb_R2(MVB)
% returns the proportion of variance explained by the (MVB) MAP estimates
% FORMAT [R2,X,P] = spm_mvb_R2(MVB)
%
% MVB - MVB structure
% R2  - proportion of variance explained
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb_R2.m 4402 2011-07-21 12:37:24Z karl $
 

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
