function [E,dE] = spm_dexpm(A,dA)
% Differentiate a matrix exponential
% FORMAT [E,dE] = spm_dexpm(A,dA)
% A  - Lie algebra
% dA - basis function to differentiate with respect to
% E  - expm(A)
% dE - (expm(A+eps*dA)-expm(A-eps*dA))/(2*eps)
%
% Note that the algorithm is a bit slow, and should perhaps be
% re-written to use eg scaling and squaring (see Moler's dubious matrix
% exponentials paper).
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dexpm.m 7408 2018-08-24 14:54:57Z john $

if nargin<2
    dA = zeros([size(A) 0]);
else
    if (size(A,1)==1 && size(A,2)==size(dA,3)) || (size(A,2)==1 && size(A,1)==size(dA,3))
        p = A(:);
        A = zeros(size(dA,1),size(dA,2));
        for m=1:size(dA,3)
            A = A + p(m)*dA(:,:,m);
        end
    end
end 
An  = A;
E   = eye(size(A))+A;

if nargout>1
    dAn = dA;
    dE  = dA;
end
for k=2:10000
    if nargout>1
        for m=1:size(dA,3)
            dAn(:,:,m) = (dAn(:,:,m)*A + An*dA(:,:,m))/k;
        end
        dE  = dE + dAn;
    end
    An  = (An*A)/k;
    E   =  E +  An;
    if sum(An(:).^2)<numel(A)*eps^2, break; end
end

