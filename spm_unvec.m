function [varargout] = spm_unvec(vX,varargin)
% unvectorises a vectorised array 
% FORMAT [X] = spm_unvec(vX,X);
% X  - numeric, cell or stucture array
% vX - spm_vec(X)
%
% i.e. X      = spm_unvec(spm_vec(X),X)
%      [X{:}] = spm_unvec(spm_vec(X{:}),X{:}) 
%                                              - (i.e. can also deal)
%
% see spm_vec
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_vec.m 184 2005-05-31 13:23:32Z karl $

% deal to multiple outputs if necessary
%--------------------------------------------------------------------------
if nargout > 1
    varargout = spm_unvec(vX,varargin);
    return
end
if length(varargin) == 1
    X = varargin{1};
else
    X = varargin;
end

% fill in structure arrays
%--------------------------------------------------------------------------
if isstruct(X)
    f = fieldnames(X);
    for i = 1:length(f)
        c          = {X.(f{i})};
        n          = length(spm_vec(c));
        c          = spm_unvec(vX(1:n),c);
        [X.(f{i})] = deal(c{:});
        vX         = vX(n + 1:end);
    end
    varargout      = {X};
    return
end

% fill in cells arrays
%--------------------------------------------------------------------------
if iscell(X)
    for i = 1:length(X(:))
        n     = length(spm_vec(X{i}));
        X{i}  = spm_unvec(vX(1:n),X{i});
        vX    = vX(n + 1:end);
    end
    varargout      = {X};
    return
end

% reshape numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X)
    if length(size(X) > 2)
        X(:) = full(vX);
    else
        X(:) = vX;
    end
else
    X     = [];
end
varargout = {X};
