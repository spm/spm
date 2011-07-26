function [B, W] = spm_robust_glm(Y, X, dim, ks)
% Apply robust GLM
% FORMAT [B, W] = spm_robust_glm(Y, X, dim, ks)
% Y      - data matrix
% X      - design matrix
% dim    - the dimension along which the function will work
% ks     - offset of the weighting function (default: 3)
%
% OUTPUT:
% B      - parameter estimates
% W      - estimated weights
%
% Implementation of:
%   Wager TD, Keller MC, Lacey SC, Jonides J.
%   Increased sensitivity in neuroimaging analyses using robust regression.
%   Neuroimage. 2005 May 15;26(1):99-113
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% James Kilner,  Vladimir Litvak
% $Id: spm_robust_glm.m 4407 2011-07-26 12:13:00Z vladimir $

if nargin < 3 || isempty(ks)
    ks = 3;
end

if nargin < 2 || isempty(dim)
    dim = 1;
end

%-Remember the original data size and size of the mean
%--------------------------------------------------------------------------
origsize       = size(Y);
borigsize      = origsize;
borigsize(dim) = size(X, 2);

%-Convert the data to repetitions x points matrix
%--------------------------------------------------------------------------
if dim > 1
    Y  = shiftdim(Y, dim-1);  
end

if length(origsize) > 2
    Y  = reshape(Y, size(Y, 1), []);
end

%-Check the design matrix and compute leverages
%--------------------------------------------------------------------------

if size(X, 1) ~= size(Y, 1)
    error('The number of rows in the design matrix should match dimension of interest.');
end

H = diag(X*inv(X'*X)*X');
H = repmat(H(:), 1, size(Y, 2));

%-Rescale the data
%--------------------------------------------------------------------------
[Y, scalefactor] = spm_cond_units(Y);

%-Actual robust GLM
%--------------------------------------------------------------------------
ores=1;
nres=10;
n=0;

YY = Y;
YY(isnan(YY)) = 0;

while max(abs(ores-nres))>sqrt(1E-8)

    ores=nres;
    n=n+1;

    if n == 1
        W = ones(size(Y));
        W(isnan(Y))  =  0;
    end
    
    B = zeros(size(X, 2), size(Y, 2));
            
    for i = 1:size(Y, 2);
        B(:, i) = inv(X'*diag(W(:, i))*X)*X'*diag(W(:, i))*YY(:, i);
    end

    if n > 200
        warning('Robust GLM could not converge. Maximal number of iterations exceeded.');
        break;
    end

    res = Y-X*B;

    mad = nanmedian(abs(res-repmat(nanmedian(res), size(res, 1), 1)));
    res = res./repmat(mad, size(res, 1), 1);
    res = res.*H;
    res = abs(res)-ks;
    res(res<0)=0;
    nres= (sum(res(~isnan(res)).^2));
    W = (abs(res)<1) .* ((1 - res.^2).^2);
    W(isnan(Y)) = 0;
    W(Y == 0)   = 0; %Assuming X is a real measurement
end

disp(['Robust GLM finished after ' num2str(n) ' iterations.']); 

%-Restore the betas and weights to the original data dimensions
%--------------------------------------------------------------------------
B = B./scalefactor;

if length(origsize) > 2   
    B  = reshape(B, circshift(borigsize, [1 -(dim-1)]));
    W  = reshape(W, circshift(origsize,  [1 -(dim-1)]));
end

if dim > 1
    B  = shiftdim(B, length(origsize)-dim+1);
    W  = shiftdim(W, length(origsize)-dim+1);
end

%-Helper function
%--------------------------------------------------------------------------
function Y = nanmedian(X)
if ~any(any(isnan(X)))
    Y = median(X);
else
    Y = zeros(1, size(X,2));
    for i = 1:size(X, 2)
        Y(i) = median(X(~isnan(X(:, i)), i));
    end
end



