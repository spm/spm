function [slice] = spm_vb_taylor_R (Y,slice)
% Get Taylor sereis approximation to posterior correlation matrices - see paper VB3
% FORMAT [slice] = spm_vb_taylor_R (Y,slice)
%
% Y     data
% slice VB-GLMAR data structure
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_taylor_R.m 1143 2008-02-07 19:33:33Z spm $

% Get mean hyperparameter values
h0=[];
if slice.p > 0
    if size(slice.ap_mean,2)==1
        % Single voxel in slice
        a=slice.ap_mean';
        a_cov=slice.a_cov{1};
    else
        a=mean(slice.ap_mean');
        a_covs=cat(3,slice.a_cov{:});
        a_cov=mean(a_covs,3);
    end
    slice.mean.a=a;
    slice.mean.a_cov=a_cov;
    h0=a';
end
slice.mean.b=mean(slice.b');

lambda=mean(slice.mean_lambda);
slice.mean.lambda=lambda;
h0=[h0;lambda];

R=spm_vb_get_R(slice,h0);
slice.mean.R=R;

delta=0.0001;
% Get first order Taylor terms about slice
% mean values of a and lambda
for i=1:length(h0),
    h=h0;
    h(i)=h(i)-delta;
    
    [R1]=spm_vb_get_R(slice,h);
    
    h=h0;
    h(i)=h(i)+delta;
    % Loop over hyperparameters
    [R2]=spm_vb_get_R(slice,h);
    
    dR(:,:,i)=(R2-R1)/(2*delta);
end
slice.mean.h0=h0;
slice.mean.dR=dR;
