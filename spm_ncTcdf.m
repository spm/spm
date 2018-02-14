function F = spm_ncTcdf(t,v,d)
% Cumulative Distribution Function (CDF) of non-central t-distribution
% FORMAT F = spm_ncTcdf(t,v,d)
% t - T-variate (Student's t has range (-Inf,Inf))
% v - degrees of freedom (v>0, non-integer d.f. accepted)
% d - non-centrality parameter
% F - CDF of non-central t-distribution with v d.f. at points t
%
% Reference:
%--------------------------------------------------------------------------
% Algorithm AS 243: Cumulative Distribution Function of the Non-Central t
% Distribution
% Russell V. Lenth
% Applied Statistics, Vol. 38, No. 1 (1989), pp. 185-189
%__________________________________________________________________________
% Copyright (C) 2009-2018 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ncTcdf.m 7258 2018-02-14 13:09:46Z guillaume $


tn = t < 0;
if any(tn)
    F = t;
    if ~all(tn), F(~tn) = spm_ncTcdf(t(~tn),v,d); end
    F(tn)  = 1 - spm_ncTcdf(-t(tn),v,-d);
    return;
end

en     = 1;
x      = t .* t ./ (t .* t + v);
lambda = d * d;
p      = 1/2 * exp(-1/2 * lambda);
q      = sqrt(2/pi) * p * d;
s      = 1/2 - p;
a      = 1/2;
b      = 1/2 * v;
rxb    = (1 - x).^b;
albeta = log(sqrt(pi)) + gammaln(b) - gammaln(a+b);
xodd   = betainc(x,a,b);
godd   = 2 * rxb .* exp(a * log(x) - albeta);
xeven  = 1 - rxb;
geven  = b * x .* rxb;
F      = p * xodd + q * xeven;

while true
    a     = a + 1;
    xodd  = xodd - godd;
    xeven = xeven - geven;
    godd  = godd .* x * (a + b - 1)/a;
    geven = geven .* x * (a + b - 1/2) / (a + 1/2);
    p     = p * lambda / (2 * en);
    q     = q * lambda / (2 * en + 1);
    s     = s - p;
    en    = en + 1;
    F     = F + p * xodd + q * xeven;
    if en > 1024 || max(2 * s * (xodd - godd)) < 1e-12
        break;
    end
end

F = F + spm_Ncdf(-d);
