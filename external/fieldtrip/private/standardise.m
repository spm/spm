function x = standardise(x,dim)

% X = STANDARDISE(X, DIM) computes the zscore of a matrix along dimension dim
% has extended functionality as compared to the stats-toolbox's zscore function

% Copyright (C) 2009, Jan-Mathijs Schoffelen
%
% $Log: standardise.m,v $
% Revision 1.1  2009/05/19 15:59:11  jansch
% first commitment into cvs
%

n      = size(x,dim);
mx     = mean(x,dim);
%sx     = std(x,0,dim);
sx     = std(x,1,dim);
repvec = ones(1,length(size(x)));
repvec(dim) = n; 
x      = (x - repmat(mx,repvec))./repmat(sx,repvec);
