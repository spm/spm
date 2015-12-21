function [Y] = spm_dot(X,x,DIM)
% Mulitimensional dot (inner) preoduct
% [Y] = spm_dot(X,x,DIM)
%
% X  - numeric array
% x  - vector or cell array of numeric vectors
%
% Y  - inner product obtined by summ the products of X and x along DIM.
%
% if DIM is not specified the last dimension of X is used.  If x is a cell
% array recursive dot products are computed (starting with the last entry 
% if (the vector) DIM is not specified
%
% See also: spm_cross
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dot.m 6652 2015-12-21 10:51:54Z karl $

% initialise X and vX
%--------------------------------------------------------------------------
if iscell(x)
    if nargin < 3
        DIM = (1:numel(x)) + numel(size(X)) - numel(x); 
    end
    for i = 1:numel(x)
      X   = spm_dot(X,x{i},DIM(i));
      DIM = DIM - 1;
    end
    Y     = X;
    return
end

% inner product
%--------------------------------------------------------------------------
if nargin < 3, DIM = numel(size(X)); end
if isvector(X), Y = X*x; return,     end
    
d        = size(X);
d(DIM)   = 1;
Y        = zeros(d);
ind      = cell(size(d));
[ind{:}] = deal(':');
sub.type = '()';

for i = 1:numel(x)
    sub.subs      = ind;
    sub.subs{DIM} = i;
    Y = Y + subsref(X,sub)*x(i);
end
Y     = squeeze(Y);

