function [R] = spm_vb_get_R (slice,h0)
% Get posterior correlation matrix for regression coefficients
% FORMAT [R] = spm_vb_get_R (slice,h0)
% 
% Y         data
% slice     data structure
% 
% R         posterior correlation matrix of regression coefficients
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_get_R.m 1143 2008-02-07 19:33:33Z spm $


lambda=h0(end);
k=slice.k;
p=slice.p;
T=slice.T;
X=slice.X;

% This is described in the section 'approximating posterior
% covariance matrices' in paper VB3
if (slice.p==0)
    A=X'*X;
else
    a=h0(1:slice.p);
    x_err         = zeros(T-p,k);
    XtilVXtil     = zeros(k);
    for t=p+1:T,
        x_err (t-p,:) = X(t,:) - a'*slice.dX(:,:,t-p);
        aa = slice.dX(:,:,t-p)'*slice.mean.a_cov;
        XtilVXtil     = XtilVXtil + aa*slice.dX(:,:,t-p);
    end
    A = x_err'*x_err + XtilVXtil; % A in eq. 58/53 paper VB1 with completed square
end

P=lambda*A+diag(slice.mean.b);
C=inv(P);
d = diag(C);
R = C./sqrt(d*d');

