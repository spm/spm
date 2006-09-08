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
% $Id: spm_orth.m 615 2006-09-08 16:16:06Z karl $

% recursive GM orthogonlisation
%--------------------------------------------------------------------------
x     = X(:,1);
for i = 2:size(X,2)
        D     = X(:,i);
        D     = D - x*(inv(x'*x)*x'*D);
        if any(D)
                x = [x D];
        end
end

% and normalisation, if requested
%--------------------------------------------------------------------------
if nargin > 1
    x = spm_en(x);
end
