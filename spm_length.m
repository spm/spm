function [n] = spm_length(X)
% length of a vectorised numeric, cell or structure array
% FORMAT [vX] = spm_length(X)
% X  - numeric, cell or stucture array[s]
% n  - length(spm_vec(X))
%
% See spm_unvec
%__________________________________________________________________________
%
% e.g.:
% spm_vec({eye(2) 3}) = [1 0 0 1 3]'
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_length.m 6110 2014-07-21 09:36:13Z karl $


%error('spm_vec.c not compiled - see Makefile')


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
