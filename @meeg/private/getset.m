function res = getset(this, parent, fieldname, ind, values)
% Generic method for getting and setting multiple fields of meeg struct
% FORMAT res = conditions(this, ind, conditionlabels)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
%$Id $

this = struct(this);

if nargin == 3 || isempty(ind)
    ind = 1:numel(getfield(struct(this), parent));
end

% Get
if nargin <= 4
    res = cell(1, length(ind));
    for i = 1:length(ind)
        res{i} = getfield(this, parent, {ind(i)}, fieldname);
    end

    if all(cellfun('isclass', res, 'double') & ~cellfun('isempty', res))
        res = [res{:}];
    end

    if iscell(res) && (length(res) == 1)
        res = res{1};
    end

    return
end

% Set
if nargin == 5
    % This might fail in some pathological cases, but not in what it's
    % supposed to be used for.
    if isnumeric(values) && (length(values) == length(ind))
        values = num2cell(values);
    end

    for i = 1:length(ind)
        if iscell(values)
            this = setfield(this, parent, {ind(i)}, fieldname, values{i});
        else
            this = setfield(this, parent, {ind(i)}, fieldname, values);
        end
    end
    res = meeg(this);
end