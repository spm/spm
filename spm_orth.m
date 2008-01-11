function x = spm_orth(X,OPT)
% recursive orthogonalization of basis functions
% FORMAT x = spm_orth(X,OPT)
%
% X   - matrix
% OPT - 'norm' for Euclidean normalisation
%
% serial orthogionalization starting with the first column
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_orth.m 1083 2008-01-11 14:12:46Z karl $

% recursive GM orthogonlisation
%--------------------------------------------------------------------------
i     = any(X);
if ~any(i)
    x = sparse(size(X,1),0);
    return
end
X     = X(:,i);
x     = X(:,1);
for i = 2:size(X,2)
        D     = X(:,i);
        D     = D - x*(inv(x'*x)*x'*D);
        if norm(D,1) > exp(-32)
                x = [x D];
        end
end

% and normalisation, if requested
%--------------------------------------------------------------------------
if nargin > 1
    x = spm_en(x);
end
