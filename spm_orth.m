function x = spm_orth(X)
% recursive orthogonalization of basis functions
% FORMAT x = spm_orth(X)
%
% serial orthogionalization starting with the first column
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$


x     = X(:,1);
for i = 2:size(X,2)
        D     = X(:,i);
        D     = D - x*(pinv(x)*D);
        if any(D)
                x = [x D];
        end
end

