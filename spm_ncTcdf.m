function F = spm_ncTcdf(x,v,d)
% Cumulative Distribution Function (CDF) of non-central t-distribution
% FORMAT F = spm_ncTcdf(x,v,d)
% x - T-variate (Student's t has range (-Inf,Inf))
% v - degrees of freedom (v>0, non-integer d.f. accepted)
% d - non-centrality parameter
% F - CDF of non-central t-distribution with v d.f. at points x
%
% Reference:
%--------------------------------------------------------------------------
% Algorithm AS 243: Cumulative Distribution Function of the Non-Central t
% Distribution
% Russell V. Lenth
% Applied Statistics, Vol. 38, No. 1 (1989), pp. 185-189
%__________________________________________________________________________
% Copyright (C) 2009-2018 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Guillaume Flandin
% $Id: spm_ncTcdf.m 7263 2018-02-21 13:30:16Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
ad = [ndims(x);ndims(v);ndims(d)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))];...
      [size(d),ones(1,rd-ad(3))];];
rs = max(as);
xa = prod(as,2)>1;
if all(xa) && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end


%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for strictly positive v. Return NaN if undefined.
md = true(size(x)) & v>0;
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Compute where defined
if ~any(md(:)), return, end

%-Compute
Q = md & x < 0;
if any(Q(:))
    if xa(1), Qx=Q; else Qx=1; end
    if xa(2), Qv=Q; else Qv=1; end
    if xa(3), Qd=Q; else Qd=1; end
    F(Q) = 1 - sub_ncTcdf(-x(Qx),v(Qv),-d(Qd));
end
Q = md & x >= 0;
if any(Q(:))
    if xa(1), Qx=Q; else Qx=1; end
    if xa(2), Qv=Q; else Qv=1; end
    if xa(3), Qd=Q; else Qd=1; end
    F(Q) = sub_ncTcdf(x(Qx),v(Qv),d(Qd));
end


%==========================================================================
function F =  sub_ncTcdf(t,v,d)
en     = 1;
x      = t .* t ./ (t .* t + v);
lambda = d .* d;
p      = 1/2 * exp(-1/2 * lambda);
q      = sqrt(2/pi) * p .* d;
s      = 1/2 - p;
a      = 1/2;
b      = 1/2 * v;
rxb    = (1 - x).^b;
albeta = log(sqrt(pi)) + gammaln(b) - gammaln(a+b);
xodd   = betainc(x,a,b);
godd   = 2 .* rxb .* exp(a * log(x) - albeta);
xeven  = 1 - rxb;
geven  = b .* x .* rxb;
F      = p .* xodd + q .* xeven;

while true
    a     = a + 1;
    xodd  = xodd - godd;
    xeven = xeven - geven;
    godd  = godd .* x .* (a + b - 1) ./ a;
    geven = geven .* x .* (a + b - 1/2) ./ (a + 1/2);
    p     = p .* lambda ./ (2 * en);
    q     = q .* lambda ./ (2 * en + 1);
    s     = s - p;
    en    = en + 1;
    F     = F + p .* xodd + q .* xeven;
    if en > 1024 || max(2 * s .* (xodd - godd)) < 1e-12
        break;
    end
end

F = F + spm_Ncdf(-d);
