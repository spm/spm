function [varargout] = spm_unvec(vX,varargin)
% Unvectorise a vectorised array - a compiled routine
% FORMAT [varargout] = spm_unvec(vX,varargin)
% varargin  - numeric, cell or structure array
% vX        - spm_vec(X)
%
% i.e. X           = spm_unvec(spm_vec(X),X)
%      [X1,X2,...] = spm_unvec(spm_vec(X1,X2,...),X1,X2,...)
%
% See spm_vec
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_unvec.m 6238 2014-10-13 09:38:23Z karl $


%error('spm_unvec.c not compiled - see Makefile')

% deal to multiple outputs if necessary
%--------------------------------------------------------------------------
if nargout > 1
    varargout = spm_unvec(vX,varargin);
    return
end
if numel(varargin) == 1
    X  = varargin{1};
else
    X  = varargin;
end

% check vX is numeric
%--------------------------------------------------------------------------
if ~(isnumeric(vX) || islogical(vX))
    vX = spm_vec(vX);
end


% reshape numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X) || islogical(X)
    if ndims(X) > 2
        X(:)  = full(vX);
    else
        X(:)  = vX;
    end
    varargout = {X};
    return
end

% fill in structure arrays
%--------------------------------------------------------------------------
if isstruct(X)
    f = fieldnames(X);
    for i = 1:numel(f)
        c          = {X.(f{i})};
        if isnumeric(c)
            n      = numel(c);
        else
            n      = spm_length(c);
        end
        c          = spm_unvec(vX(1:n),c);
        [X.(f{i})] = deal(c{:});
        vX         = vX(n + 1:end);
    end
    varargout      = {X};
    return
end

% fill in cell arrays
%--------------------------------------------------------------------------
if iscell(X)
    for i = 1:numel(X)
        if isnumeric(X{i})
            n      = numel(X{i});
        else
            n      = spm_length(X{i});
        end
        X{i}  = spm_unvec(vX(1:n),X{i});
        vX    = vX(n + 1:end);
    end
    varargout = {X};
    return
end

% else
%--------------------------------------------------------------------------
X         = [];
varargout = {X};
