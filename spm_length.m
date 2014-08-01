function [n] = spm_length(X)
% Length of a vectorised numeric, cell or structure array
% FORMAT [n] = spm_length(X)
% X    - numeric, cell or stucture array[s]
% n    - length(spm_vec(X))
%
% See spm_vec, spm_unvec
%__________________________________________________________________________
%
% e.g.:
% spm_length({eye(2) 3}) = 5
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_length.m 6130 2014-08-01 17:41:18Z guillaume $


% vectorise numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X) 
    n = numel(X);

% vectorise logical arrays
%--------------------------------------------------------------------------
elseif islogical(X)
    n = numel(X);

% vectorise structure into cell arrays
%--------------------------------------------------------------------------
elseif isstruct(X)
    n     = 0;
    f     = fieldnames(X);
    for i = 1:numel(f)
        n = n + spm_length(X.(f{i}));
    end

% vectorise cells into numerical arrays
%--------------------------------------------------------------------------
elseif iscell(X)
    n     = 0;
    for i = 1:numel(X)
        n = n + spm_length(X{i});
    end
    
else
    n = 0;
end
