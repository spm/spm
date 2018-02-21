function f = spm_ncFpdf(x,df,d)
% Probability Density Function (PDF) of non-central F-distribution
% FORMAT f = spm_ncFpdf(x,df,d)
% x  - F-variate (F has range [0,Inf) )
% df - degrees of freedom, df = [v,w] with v>0 and w>0
% d  - non-centrality parameter
% f  - PDF of non-central F-distribution with [v,w] degrees of freedom (df)
%      and non-centrality parameter d at points x
%__________________________________________________________________________
%
% spm_ncFpdf implements the Probability Density Function of non-central 
% F-distributions.
%
% References:
% https://en.wikipedia.org/wiki/Noncentral_F-distribution
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_ncFpdf.m 7263 2018-02-21 13:30:16Z guillaume $


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
f = zeros(rs);

%-Only defined for v>0, w>0 and d>=0. Return NaN if undefined.
md = true(size(f)) & v>0 & w>0 & d>=0;
if any(~md(:))
    f(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Compute where defined
if ~any(md(:)), return, end

Q = md & x>=0; % 0 if x < 0
if any(Q(:))
    if xa(1), Qx=Q; else Qx=1; end
    if xa(2), Qv=Q; else Qv=1; end
    if xa(3), Qw=Q; else Qw=1; end
    if xa(4), Qd=Q; else Qd=1; end
    
    a = (v(Qv) ./ w(Qw)).^(v(Qv)/2) ./ ...
        beta(w(Qw)/2,v(Qv)/2) .* ...
        (w(Qw) ./ (w(Qw) + v(Qv) .* x(Qx))).^((v(Qv)+w(Qw))/2) .* ...
        x(Qx).^(v(Qv)/2-1);
    f(Q) = a;
    for i=1:1024
        % Use recursion B(x,y+k) = B(x,y+k-1) * (y+k-1) / (x+y+k-1)
        a = a .* (d(Qd)/2) ./ ...
            ((v(Qv)/2+i-1) ./ (v(Qv)/2+w(Qw)/2+i-1) * i) .* ...
            (v(Qv) ./ w(Qw)) .* ...
            (w(Qw) ./ (w(Qw) + v(Qv) .* x(Qx))) .* ...
            x(Qx);
        f(Q) = f(Q) + a;
        if max(a) < 1e-12, break; end
    end
    f(Q) = exp(-d(Qd)/2) .* f(Q);
end
