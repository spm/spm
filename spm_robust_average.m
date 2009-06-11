function [Y,W] = spm_eeg_robust_average(X, dim, ks, h)
% Apply robust averaging routine to X sets
% FORMAT [Y,W] = spm_eeg_robust_averaget(X, dim, ks)
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
% $Id: spm_robust_average.m 3196 2009-06-11 12:54:47Z vladimir $

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
Y=zeros(1,size(X,1));
while max(abs(ores-nres))>sqrt(1E-8)

    ores=nres;
    n=n+1;

    if n==1
        Y = median(X);
    else
        Y = sum(W.*X)./sum(W);
    end

    if sum(isnan(Y))>0
        error('NaNs appeared in the result');
    end

    if n > 100
        warning('Robust averaging could not converge. Maximal number of iterations exceeded.');
    end

    res = X-repmat(Y, size(X, 1), 1);

    mad = median(abs(res-repmat(median(res), size(res, 1), 1)));
    res = res./repmat(mad, size(res, 1), 1);
    res = res.*h;
    res = abs(res)-ks;
    res(res<0)=0;
    nres= (sum(res(:).^2));
    W = (abs(res)<1) .* ((1 - res.^2).^2);
end

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


