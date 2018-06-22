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
% $Id: spm_ncFcdf.m 7354 2018-06-22 10:44:22Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
%-Unpack degrees of freedom v & w from single df parameter (df)
if numel(df) == 2
    v = df(1);
    w = df(2);
elseif size(df,2) == 2
    v = df(:,1);
    w = df(:,2);
else
    error('Cannot unpack degrees of freedom.');
end

%-Check argument sizes
ad = [ndims(x);ndims(v);ndims(w);ndims(d)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))];...
      [size(w),ones(1,rd-ad(3))];...
      [size(d),ones(1,rd-ad(4))];];
rs = max(as);
xa = prod(as,2)>1;
if all(xa) && any(any(diff(as(xa,:)),1))
    error('Non-scalar arguments must match in size.');
end


%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for v>0, w>0 and d>=0. Return NaN if undefined.
md = true(size(F)) & v>0 & w>0 & d>=0;
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Compute where defined
if ~any(md(:)), return, end
Q = md;
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end
if xa(4), Qd=Q; else Qd=1; end

a = exp(-d(Qd)/2);
x = v(Qv) .* x(Qx) ./ (w(Qw) + v(Qv) .* x(Qx));
F(Q) = a .* betainc(x(Qx), v(Qv)/2, w(Qw)/2);
for i=1:1024
    a = a .* d(Qd) / (2 * i);
    % could use recursion with:
    % Ix(a+1,b)=Ix(a,b)-G(a+b)/(G(a+1)*G(b))*x^a*(1-x)^b and G(a+1)=a*G(a)
    e = a .* betainc(x(Qx), v(Qv)/2+i, w(Qw)/2);
    if max(e) < 1e-12, break; end
    F = F + e;
end
