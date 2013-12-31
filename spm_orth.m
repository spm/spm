function X = spm_orth(X,OPT)
% Recursive Gram-Schmidt orthogonalisation of basis functions
% FORMAT X = spm_orth(X,OPT)
%
% X   - matrix
% OPT - 'norm' for Euclidean normalisation
%     - 'pad'  for zero padding of null space [default]
%
% Serial orthogonalisation starting with the first column
%
% Reference:
% Golub, Gene H. & Van Loan, Charles F. (1996), Matrix Computations (3rd
% ed.), Johns Hopkins, ISBN 978-0-8018-5414-9.
%__________________________________________________________________________
% Copyright (C) 2002-2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_orth.m 5821 2013-12-31 14:26:41Z karl $
 

%-Default
%--------------------------------------------------------------------------
try
    OPT;
catch
    OPT = 'pad';
end
 
%-Recursive Gram-Schmidt orthogonalisation
%--------------------------------------------------------------------------
sw    = warning('off','all');
[n,m] = size(X);
X     = X(:, any(X));
rankX = rank(full(X));
try
    x     = X(:,1);
    j     = 1;
    for i = 2:size(X, 2)
        D = X(:,i);
        D = D - x*(pinv(x)*D);
        if norm(D,1) > exp(-32)
            x          = [x D];
            j(end + 1) = i;
        end
        if numel(j) == rankX, break, end
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
        X      = zeros(n,m);
        X(:,j) = x;
    case{'norm'}
        X      = spm_en(x);
    otherwise
        X      = x;
end
