function F = spm_ncFcdf(x,df,d)
% Cumulative Distribution Function (CDF) of non-central F-distribution
% FORMAT f = spm_ncFcdf(x,df,d)
% x  - F-variate (F has range [0,Inf) )
% df - degrees of freedom, df = [v,w] with v>0 and w>0
% d  - non-centrality parameter
% F  - CDF of non-central F-distribution with [v,w] d.f. at points x
%
% Reference:
% https://en.wikipedia.org/wiki/Noncentral_F-distribution
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_ncFcdf.m 7263 2018-02-21 13:30:16Z guillaume $


if any(d(:) < 0) || any(df(:) <= 0)
    F = NaN(size(x));
    warning('Returning NaN for out of range arguments');
    return;
end

a = exp(-d/2);
x = df(1) * x ./ (df(2) + df(1) * x);
F = a .* betainc(x, df(1)/2, df(2)/2);
for i=1:1024
    a = a .* d / (2 * i);
    % could use recursion with:
    % Ix(a+1,b)=Ix(a,b)-G(a+b)/(G(a+1)*G(b))*x^a*(1-x)^b and G(a+1)=a*G(a)
    e = a .* betainc(x, df(1)/2+i, df(2)/2);
    if max(e) < 1e-12, break; end
    F = F + e;
end
