function X = spm_orth(X,OPT)
% recursive orthogonalisation of basis functions
% FORMAT X = spm_orth(X,OPT)
%
% X   - matrix
% OPT - 'norm' for Euclidean normalisation
%     - 'pad'  for zero padding of null space [default]
%
% serial orthogonalisation starting with the first column
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_orth.m 1172 2008-02-27 20:14:47Z karl $
 
% default
%--------------------------------------------------------------------------
try
    OPT;
catch
    OPT = 'pad';
end
 
% recursive GM orthogonalisation
%--------------------------------------------------------------------------
sw = warning('off','all');
[n m] = size(X);
i     = find(any(X));
X     = X(:,i);
try
    x     = X(:,1);
    j     = 1;
    for i = 2:size(X,2)
        D = X(:,i);
        D = D - x*(inv(x'*x)*x'*D);
        if norm(D,1) > exp(-32)
            x          = [x D];
            j(end + 1) = i;
        end
    end
catch
    x     = zeros(n,0);
    j     = [];
end
warning(sw);
 
% and normalisation, if requested
%--------------------------------------------------------------------------
switch OPT
    case{'pad'}
        X      = sparse(n,m);
        X(:,j) = x;
    otherwise
        X      = spm_en(x);
end

