function res = getset(this, parent, fieldname, ind, values)
% Generic method for getting and setting multiple fields of meeg struct
% FORMAT res = getset(this, parent, fieldname, ind, values)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


this = struct(this);

if nargin == 3 || ~isnumeric(ind)
    try
        ind = 1:numel(getfield(this, parent));
    catch
        res = [];
        return;
    end
end

% Get
if nargin <= 4
    try
        res = {this.(parent)(ind).(fieldname)};
    catch
        res = cell(1, length(ind));
        for i = 1:length(ind)
            res{i} = getfield(this, parent, {ind(i)}, fieldname);
        end
    end

    if isempty(res) || (all(cellfun('isclass', res, 'double') & (cellfun(@numel, res) == 1)))
        res = [res{:}];
    end

    if iscell(res) && (numel(res) == 1) && (numel(getfield(this, parent)) == 1) &&...
            strcmp(parent, 'trials') && strcmp(this.type, 'continuous')
        res = res{1};
    end
    
    return
end

% Set
if nargin == 5          
    % This might fail in some pathological cases, but not in what it's
    % supposed to be used for.
    if (isnumeric(values) || islogical(values)) && (length(values) == length(ind))
        values = num2cell(values);
    end
    
    if iscell(values) && ~(numel(values)==1 || numel(values) == length(ind))
        error('Illegal assignment: cannot match values and indices.');
    end        
        
    if iscell(values)
        try 
            [this.(parent)(ind).(fieldname)] = values{:};  
        catch 
            this = setfield(this, parent, {ind(i)}, fieldname, values{i});
        end
    else 
        try 
            [this.(parent)(ind).(fieldname)] = deal(values);  
        catch 
            this = setfield(this, parent, {ind(i)}, fieldname, values);
        end
    end
    
    % getset is sometimes used on subfields of meeg then checkmeeg should
    % not be used
    if all(isfield(this, {'type', 'Nsamples', 'Fsample', 'timeOnset'}))
        res = meeg(this);
    else
        res = this;
    end
end
