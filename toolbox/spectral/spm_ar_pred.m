function [y_pred,y,r2] = spm_ar_pred (Z,ar)
% Make predictions from Bayesian autoregressive models
% FORMAT [y_pred,y,r2] = spm_ar_pred (Z,ar)
%
% Z             [N x 1] univariate time series 
% ar            data structure - see spm_ar.m
%
% y_pred        (one-step ahead) predictions 
% y             the values we are 'predicting'
% r2            proportion of variance explained
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_ar_pred.m 1143 2008-02-07 19:33:33Z spm $

if nargin < 2, 
   disp('spm_ar_pred.m needs at least two arguments'); 
   return
end

p=ar.p;

Z=Z(:);
y=Z(p+1:end);
for i=1:p,
    x(:,i)=Z(p-i+1:end-i);
end

y_pred = -x*ar.a_mean;

vy=std(y)^2;
ey=std(y-y_pred)^2;
r2=(vy-ey)/vy; 