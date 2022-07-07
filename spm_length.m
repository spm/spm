function [n] = spm_length(X)
% Length of a vectorised numeric, cell or structure array
% FORMAT [n] = spm_length(X)
% X    - numeric, cell or structure array[s]
% n    - length(spm_vec(X))
%
% See spm_vec, spm_unvec
%__________________________________________________________________________
%
% e.g.:
% spm_length({eye(2) 3}) = 5
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2014-2022 Wellcome Centre for Human Neuroimaging


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
        for j = 1:numel(X)
            n = n + spm_length(X(j).(f{i}));
        end
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
