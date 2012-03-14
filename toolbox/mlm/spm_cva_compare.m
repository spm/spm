function [CVA] = spm_cva_compare (Y,X,c)
% Model comparison for probabilistic CVA
% FORMAT [CVA] = spm_cva_compare (Y,X,c)
%
% Y  [N x d1] data matrix
% X  [N x d2] design matrix
% c  Contrast vector (if specified)
%
% CVA has fields:
%
% .order        number of canonical vectors (latent space dimension)
% .bic          BIC for each order
% .aic          AIC for each order
%
% and 
%
% .U1,.U2      Canonical vectors
% .W1,.W2      Factor matrices
%
% for the highest order model.
%
% See spm_cva_prob.m for more details
%___________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_cva_compare.m 4687 2012-03-14 18:15:49Z will $

if nargin > 2
    % get null-space of contrast
    X0  = X - X*c*pinv(c);
    X   = full(X*c);
    X0  = spm_svd(X0);
    
    Y     = Y - X0*(X0'*Y);
    X     = X - X0*(X0'*X);
end

[N1,p1]=size(Y);
[N2,p2]=size(X);
if ~(N1==N2)
    disp('X and Y are of incompatible size');
    return
end

m=min([p1,p2]);
for i=1:m+1,
    order(i)=i-1;
    CVA = spm_cva_prob (Y',X',order(i));
    bic(i)=CVA.bic;
    aic(i)=CVA.aic;
    L(i)=CVA.L;
end

CVA.order=order;
CVA.bic=bic;
CVA.aic=aic;
CVA.L=L;
