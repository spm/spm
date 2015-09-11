function [p] = spm_mvtpdf (x,mu,Lambda,v)
% PDF of multivariate T-distribution
% FORMAT [p] = spm_mvtpdf (x,mu,Lambda,v)
%
% x      - ordinates [d x N]
% mu     - mean [d x 1]
% Lambda - precision matrix [d x d]
% v      - degrees of freedom
%
% p      - probability density
%
% See J. Bernardo and A. Smith (2000) 
% Bayesian Theory, Wiley (page 435)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mvtpdf.m 6548 2015-09-11 12:39:47Z will $

d=length(mu);
lnz=gammaln(0.5*v)+0.5*d*log(v)+0.5*d*log(pi);
lnz=lnz-0.5*spm_logdet(Lambda)-gammaln(0.5*(v+d));
z=exp(lnz);

N=size(x,2);
for n=1:N,
    e=x(:,n)-mu;
    p(n)=(1+(1/v)*e'*Lambda*e).^(-0.5*(v+d));
end
p=p/z;
