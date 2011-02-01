function [p] = spm_ncTcdf(x,v,d)
% Cumulative Distribution Function (CDF) of non-central t-distribution
% FORMAT [p] = spm_ncTcdf(x,v,d)
% x -  T-variate (Student's t has range (-Inf,Inf)
% v - degrees of freedom (v>0, non-integer d.f. accepted)
% F - CDF of non-central t-distribution with v d.f. at points x
%
% see:
% Algorithm AS 243: Cumulative Distribution Function of the Non-Central t
% Distribution
% Russell V. Lenth
% Applied Statistics, Vol. 38, No. 1 (1989), pp. 185-189
%__________________________________________________________________________
% Copyright (C) 2009-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ncTcdf.m 4182 2011-02-01 12:29:09Z guillaume $


% CDF
%--------------------------------------------------------------------------
t     = x;
p     = x;
ip    = find(x >= 0);
in    = find(x <  0);


% positive support
%--------------------------------------------------------------------------
x     = t(ip);
x     = x.^2./(x.^2 + v);
P     = 0;
for j = 0:128
    Ix   = betainc(x,j + 1/2,v/2);
    Jx   = betainc(x,j + 1,v/2);
    pj   =   (1/2)*exp(-d^2/2)*(d^2/2)^j/factorial(j);
    qj   = d*(1/2)*exp(-d^2/2)*(d^2/2)^j/(sqrt(2)*gamma(j + 3/2));
    P    = P + pj*Ix + qj*Jx;
end
p(ip)    = spm_Ncdf(-d) + P;


% and negative support
%--------------------------------------------------------------------------
d     = -d;
x     = t(in);
x     = x.^2./(x.^2 + v);
P     = 0;
for j = 0:32
    Ix   = betainc(x,j + 1/2,v/2);
    Jx   = betainc(x,j + 1,v/2);
    pj   =   (1/2)*exp(-d^2/2)*(d^2/2)^j/factorial(j);
    qj   = d*(1/2)*exp(-d^2/2)*(d^2/2)^j/(sqrt(2)*gamma(j + 3/2));
    P    = P + pj*Ix + qj*Jx;
end
p(in)    = 1 - spm_Ncdf(-d) - P;
