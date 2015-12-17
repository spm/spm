function [Y,W] = spm_robust_average(X, dim, ks)
% Apply robust averaging routine to X sets
% FORMAT [Y,W] = spm_robust_averaget(X, dim, ks)
% X      - data matrix to be averaged
% dim    - the dimension along which the function will work
% ks     - offset of the weighting function (default: 3)
%
% W      - estimated weights
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% James Kilner
% $Id: spm_robust_average.m 6650 2015-12-17 14:48:43Z vladimir $

if nargin < 3 || isempty(ks)
    ks = 3;
end

if nargin < 2 || isempty(dim)
    dim = 1;
end

%-Remember the original data size and size of the mean
%--------------------------------------------------------------------------
origsize       = size(X);
morigsize      = origsize;
morigsize(dim) = 1;

if length(origsize)<dim || origsize(dim) == 1
    warning('There is only one replication in the data. Robust averaging cannot be done.');
    Y = X;
    W = ones(size(X));
    return;
end

%-Convert the data to repetitions x points matrix
%--------------------------------------------------------------------------
if dim > 1
    X  = shiftdim(X, dim-1);  
end

if length(origsize) > 2
    X  = reshape(X, size(X, 1), []);
end

%-Replace Inf with NaN
%--------------------------------------------------------------------------
X(~isfinite(X)) = NaN;

%-Rescale the data
%--------------------------------------------------------------------------
[X, scalefactor] = spm_cond_units(X);

%-Actual robust averaging
%--------------------------------------------------------------------------
ores=1;
nres=10;
n=0;

W = zeros(size(X));

while max(abs(ores-nres))>sqrt(1E-8)

    ores=nres;
    n=n+1;

    if n==1
            Y = nanmedian(X); 
    else
        XX = X;
        XX(isnan(XX)) = 0;
        Y = sum(W.*XX)./sum(W);
    end

    if n > 200
        warning('Robust averaging could not converge. Maximal number of iterations exceeded.');
        break;
    end

    res = X-repmat(Y, size(X, 1), 1);
    
    mad = nanmedian(abs(res-repmat(nanmedian(res), size(res, 1), 1)));
    
    ind1 = find(mad==0);
    ind2 = find(mad~=0);
    
    res1       = res(:, ind1);
    res1(isnan(res1)) = 1;
    W(:, ind1) = ~res1;
    
    if ~isempty(ind2)
        res = res(:, ind2);
        mad = mad(ind2);
        
        res = res./repmat(mad, size(res, 1), 1);
        res = abs(res)-ks;
        res(res<0) = 0;
        nres = (sum(res(~isnan(res)).^2));
        W(:, ind2)  = (abs(res)<1) .* ((1 - res.^2).^2);
        W(W == 0) = eps; % This is to prevent appearance of NaNs when normalizing
        W(isnan(X)) = 0;
        W(X == 0 & ~repmat(all(X==0), size(X, 1), 1)) = 0; %Assuming X is a real measurement        
    end      
end

disp(['Robust averaging finished after ' num2str(n) ' iterations.']); 

%-Restore the average and weights to the original data dimensions
%--------------------------------------------------------------------------
Y = Y./scalefactor;

if length(origsize) > 2   
    Y  = reshape(Y, circshift(morigsize, [1 -(dim-1)]));
    W  = reshape(W, circshift(origsize,  [1 -(dim-1)]));
end

if dim > 1
    Y  = shiftdim(Y, length(origsize)-dim+1);
    W  = shiftdim(W, length(origsize)-dim+1);
    
    % This is helpful when there are singleton dimensions
    Y  = reshape(Y, morigsize);
    W  = reshape(W, origsize);
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



