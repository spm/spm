function f = spm_ncTpdf(x,v,d)
% Probability Density Function (PDF) of non-central t-distribution
% FORMAT f = spm_ncTpdf(x,v,d)
% x - t-ordinates
% v - degrees of freedom (v>0, non-integer d.f. accepted)
% d - non-centrality parameter
% f - PDF of non-central t-distribution with v degrees of freedom (df) and
%     non-centrality parameter d at points x
%__________________________________________________________________________
%
% spm_ncTpdf implements the Probability Density Function of non-central 
% t-distributions.
%
% References:
% https://en.wikipedia.org/wiki/Noncentral_t-distribution
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_ncTpdf.m 7263 2018-02-21 13:30:16Z guillaume $


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
    error('Non-scalar arguments must match in size.');
end


%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
f = zeros(rs);

%-Only defined for v>0. Return NaN if undefined.
md = true(size(f)) & v>0;
if any(~md(:))
    f(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Compute where defined
if ~any(md(:)), return, end

%-Compute
Q = md & x == 0;
if any(Q(:))
    if xa(2), Qv=Q; else Qv=1; end
    if xa(3), Qd=Q; else Qd=1; end
    f(Q) = exp(-d(Qd).^2/2) ./ beta(v(Qv)/2,1/2) ./ sqrt(v(Qv));
end

Q = md & x ~= 0;
if any(Q(:))
    if xa(1), Qx=Q; else Qx=1; end
    if xa(2), Qv=Q; else Qv=1; end
    if xa(3), Qd=Q; else Qd=1; end
    f(Q) = v(Qv) ./ x(Qx) .* ...
        (spm_ncTcdf(x(Qx).*sqrt(1+2./v(Qv)),v(Qv)+2,d(Qd)) - spm_ncTcdf(x(Qx),v(Qv),d(Qd)));
end
