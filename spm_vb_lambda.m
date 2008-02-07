function [slice] = spm_vb_lambda (Y,slice)
% Variational Bayes for GLM-AR models - Update lambda
% FORMAT [slice] = spm_vb_lambda (Y,slice)
%
% Y             [T x N] time series 
% slice         data structure containing the following fields:
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_lambda.m 1143 2008-02-07 19:33:33Z spm $

if slice.verbose
    disp('Updating lambda');
end

for n=1:slice.N,
    if slice.p > 0
        % Equation 77 in paper VB1
        Gn = spm_vb_get_Gn (Y,slice,n);
    else
        block_n           = [(n-1)*slice.k+1:n*slice.k];
        en=(Y(:,n)-slice.X*slice.w_mean(block_n,1));
        Gn=trace(slice.w_cov{n}*slice.XTX)+en'*en;
    end
    % Equation 75 in paper VB1
    slice.b_lambda(n,1)    = 1./(Gn./2 + 1./slice.b_lambda_prior(n));
    slice.mean_lambda(n,1) = slice.c_lambda(n)*slice.b_lambda(n);
end

