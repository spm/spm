function [X] = spm_vecfun(X,fun)
% Apply a function to the numeric elements of a cell or structure array
% FORMAT [X] = spm_vecfun(X,fun)
% X   - numeric, cell or stucture array
% fun - function handle
%__________________________________________________________________________
%
% e.g., pE = spm_vecfun(pE,@log)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


% vectorise numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X)
    X = fun(X);

% vectorise structure into cell arrays
%--------------------------------------------------------------------------
elseif isstruct(X)
    f     = fieldnames(X);
    for i = 1:numel(f)
        X.(f{i}) = spm_vecfun(X.(f{i}),fun);
    end

% vectorise cells into numerical arrays
%--------------------------------------------------------------------------
elseif iscell(X)
    for i = 1:numel(X)
        X{i} = spm_vecfun(X{i},fun);
    end
end
