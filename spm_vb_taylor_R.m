function [slice] = spm_vb_taylor_R (Y,slice)
% Get Taylor sereis approximation to posterior correlation matrices - see paper VB3
% FORMAT [slice] = spm_vb_taylor_R (Y,slice)
%
% Y     data
% slice VB-GLMAR data structure
%
% @(#)spm_vb_taylor_R.m	1.1 Will Penny 04/08/04

% Get mean hyperparameter values
h0=[];
if slice.p > 0
    a=mean(slice.ap_mean');
    a_covs=cat(3,slice.a_cov{:});
    a_cov=mean(a_covs,3);
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
