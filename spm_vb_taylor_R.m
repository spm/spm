function [slice] = spm_vb_taylor_R(Y,slice)
% Get Taylor series approximation to posterior correlation matrices
% FORMAT [slice] = spm_vb_taylor_R(Y,slice)
%
% Y        - data
% slice    - VB-GLMAR data structure
%
% See paper VB3.
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_taylor_R.m 6079 2014-06-30 18:25:37Z spm $

% Get mean hyperparameter values
h0 = [];
if slice.p > 0
    if size(slice.ap_mean,2)==1
        % Single voxel in slice
        a      = slice.ap_mean';
        a_cov  = slice.a_cov{1};
    else
        a      = mean(slice.ap_mean');
        a_covs = cat(3,slice.a_cov{:});
        a_cov  = mean(a_covs,3);
    end
    slice.mean.a     = a;
    slice.mean.a_cov = a_cov;
    h0 = a';
end
slice.mean.b = mean(slice.b');

lambda = mean(slice.mean_lambda);
slice.mean.lambda = lambda;
h0 = [h0;lambda];

R = spm_vb_get_R(slice,h0);
slice.mean.R = R;

delta = 0.0001;
% Get first order Taylor terms about slice
% mean values of a and lambda
for i=1:length(h0)
    h    = h0;
    h(i) = h(i)-delta;
    
    [R1] = spm_vb_get_R(slice,h);
    
    h    = h0;
    h(i) = h(i) + delta;
    % Loop over hyperparameters
    R2   = spm_vb_get_R(slice,h);
    
    dR(:,:,i) = (R2-R1) / (2*delta);
end
slice.mean.h0 = h0;
slice.mean.dR = dR;
