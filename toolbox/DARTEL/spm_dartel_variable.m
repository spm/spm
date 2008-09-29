function [U,ll,f1] = spm_dartel_variable(U,f,g,param)
% Iteration of variable velocity DARTEL
% FORMAT [U,f1] = spm_dartel_variable(U,f,g,param)
%   U     - Flow fields. Single precision d1*d2*d3*3*k1
%   f     - Individual's data. Single precision d1*d2*d3*n
%   g     - Template. Single precision d1*d2*d3*n
%   param - Parameters (similar to those for dartel3).
%
%   ll    - Various metrics
%   f1    - Part of sufficient statistics for re-building template
%
% The idea is to approximate Beg's LDDMM algorithm by composing a
% series of DARTEL steps. This is a subroutine of that procedure.
% See also spm_dartel_integrate.m.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_variable.m 2231 2008-09-29 15:01:05Z guillaume $
 
dm         = size(U);
k1         = dm(4);
dm         = dm(1:3);
param(2:4) = param(2:4)*k1;
param(8)   = ceil(log2(2^param(8)/k1));

% Generate a sequence of warped templates
g1    = cell(1,k1);
g1{1} = g;
if k1==1,
    for i=2:k1, g1{i} = g; end
else
    g = log(max(cat(4,g,1-sum(g,4)),1e-12));
    for i = 1:(k1-1),
        v       = squeeze(single(U(:,:,:,k1-i+1,:)));
        phi0    = dartel3('Exp',v,[param(8) 1 0]); drawnow
        clear v
        if i==1,
            phi = phi0;
        else
            phi = dartel3('comp',phi0,phi); drawnow
        end
        clear phi0

        tmp = exp(dartel3('samp',g,phi)); drawnow
        sg  = sum(tmp,4);
        g1{i+1} = zeros([dm(1:3),size(tmp,4)-1],'single');
        for n=1:(size(tmp,4)-1),
            g1{i+1}(:,:,:,n) = tmp(:,:,:,n)./sg;
        end
        clear sg tmp
    end;
    clear phi
end

er2  = 0;
len2 = 0;
len1 = 0;
f1   = f;
for i = k1:-1:1,
    v         = squeeze(single(U(:,:,:,k1-i+1,:)));
    [v,ll]    = dartel3(v,f1,g1{i},param); drawnow
   %fprintf('%g\t%g\t%g\n', ll);

    er2  = er2 +ll(1)/k1;
    len2 = len2+ll(2)/k1;
    len1 = len1+sqrt(ll(2)/k1);

    if i>1 || nargout>=2,
        phi0 = dartel3('Exp',v,[param(8) 1 1]); drawnow
        U(:,:,:,k1-i+1,:) = reshape(v,[dm 1 3]);
        clear v
        if i==k1,
            phi = phi0;
        else
            phi = dartel3('comp',phi,phi0); drawnow
        end
        clear phi0 
        [f1,dJ] = dartel3('pushc',f,phi); drawnow
        if i>1,
            for n=1:size(f1,4),
                f1(:,:,:,n) = f1(:,:,:,n)./(dJ+eps);
            end
        else
            f1 = cat(4,f1,dJ);
        end
    end
end
%fprintf('*** %g %g %g %g ***\n', er2, len2, er2+len2, len1);
ll = [er2 len2 len1];
