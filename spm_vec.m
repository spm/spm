function [x] = spm_vec(varargin)
% FORMAT [x] = spm_vec(x);
% vectorises an array or structure of hierarchically embedded matrices
% x - cell array
% x - vec(x)
%____________________________________________________________________________
%
% e.g.:
% spm_vec({eye(2) 3}) = [1 0 0 1 3]'
%____________________________________________________________________________
% %W% Karl Friston %E%

% de-reference input cell if possible
%----------------------------------------------------------------------------
x = varargin;
if length(x) == 1
    x = x{1};
end

% vectorise structure into cell arrays
%----------------------------------------------------------------------------
if isstruct(x)
    f = fieldnames(x);
    y = f(:);
    for i = 1:length(f)
         y{i} = spm_vec(getfield(x,f{i}));
    end
    x = y;
end

% vectorise cells into numberical arrays
%----------------------------------------------------------------------------
if iscell(x)
    x     = x(:);
    for i = 1:length(x)
         x{i} = spm_vec(x{i});
    end
    x         = spm_cat(x);
end

% vectorise numerical arrays
%----------------------------------------------------------------------------
if isnumeric(x)
    x = x(:);
else
    x = [];
end



