function [Y,W] = spm_robust_average(X, dim, ks, h)
% Apply robust averaging routine to X sets
% FORMAT [Y,W] = spm_robust_averaget(X, dim, ks)
% X      - data matrix to be averaged
% dim    - the dimension along which the function will work
% ks     - offset of the weighting function (default: 3)
% h      - (optional) matrix of 'leverages' indicating how
%          important each data point is. This is for computations more
%          complicated than averaging e.g. regression
%
% W      - estimated weights
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% James Kilner
% $Id: spm_robust_average.m 3205 2009-06-16 10:15:00Z vladimir $

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

%-Convert the data to repetitions x points matrix
%--------------------------------------------------------------------------
if dim > 1
    X  = shiftdim(X, dim-1);  
end

if length(origsize) > 2
    X  = reshape(X, size(X, 1), []);
end

%-Check the leverages and reshape if necessary
%--------------------------------------------------------------------------
if nargin > 3
    if numel(h)>1
        if morigsize(end) == 1
            if ~isequal(morigsize(1:(end-1)), size(h))
                error('The leverages do not match the data.');
            end
            
            if dim > 2
                h  = shiftdim(h, dim-2);
            end
        else
            if ~isequal(morigsize, size(h))
                error('The leverages do not match the data.');
            end 
            
            if dim > 1
                h  = shiftdim(h, dim-1);
            end
        end
        if length(origsize) > 2
            h  = reshape(h, 1, []);
        end
    end
else
    h = 1;
end

%-Rescale the data
%--------------------------------------------------------------------------
[X, scalefactor] = spm_cond_units(X);


%-Actual robust averaging
%--------------------------------------------------------------------------
ores=1;
nres=10;
n=0;

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
    res = res./repmat(mad, size(res, 1), 1);
    res = res.*h;
    res = abs(res)-ks;
    res(res<0)=0;
    nres= (sum(res(~isnan(res)).^2));
    W = (abs(res)<1) .* ((1 - res.^2).^2);
    W(isnan(X)) = 0;
    W(X == 0)   = 0; %Assuming X is a real measurement
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



