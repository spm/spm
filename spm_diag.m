function D = spm_diag(varargin)
% Diagonal matrices and diagonals of a matrix
%
% SPM_DIAG generalises the function "diag" to also work with cell arrays.
% See DIAG's help for syntax.
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_diag.m 4310 2011-04-18 16:07:35Z guillaume $

try
    X = varargin{1};
catch
    error('Not enough input arguments.');
end

% use built-in diag for most data types
%--------------------------------------------------------------------------
if ~iscell(X)

    D = diag(varargin{:});

% or use the following for cell arrays
%--------------------------------------------------------------------------
else

    try, K = varargin{2}; catch, K = 0; end
    [m,n]  = size(X);

    % return the cell array whose K-th diagonal is X
    %----------------------------------------------------------------------
    if any([m n] == 1)
        D  = cell(max(m,n) + abs(K));
        D(logical(diag(ones(1,max(m,n)),K))) = X;

    % return the K-th diagonal of X
    %----------------------------------------------------------------------
    else
        D  = X(logical(diag(ones(1,max(m,n)-abs(K)),K)));
    end

end
