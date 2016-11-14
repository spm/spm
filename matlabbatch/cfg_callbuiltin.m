function varargout = cfg_callbuiltin(varargin)

if exist('OCTAVE_VERSION', 'builtin') && ismember(varargin{1}, {'subsref','subsasgn'})
    % Fix for inheritance in GNU Octave
    % Many thanks to Stanislaw Adaszewski (http://algoholic.eu)
    obj  = varargin{2};
    subs = varargin{3};
    varargin{3} = struct('type', {}, 'subs', {});
    for i=1:numel(subs)
        switch subs(i).type
            case '.'
                sub = findfield(obj, subs(i).subs);
            otherwise
                sub = subs(i);
        end
        try, obj = builtin('subsref', obj, sub); end
        varargin{3} = [varargin{3} sub];
    end
end

if nargout == 0
    try
        varargout{1} = builtin(varargin{:});
    catch
        builtin(varargin{:});
    end
else
    varargout = cell(1,nargout);
    [varargout{:}] = builtin(varargin{:});
end


function S = findfield(obj, field)
S = substruct('.', field);
if ~isobject(obj), return; end
fields = builtin('fieldnames', obj);
if ismember(field, fields), return; end
for i=1:numel(fields)
    if strncmp(fields{i}, 'cfg_', 4) % to check fields{i} is a classname
        sub1 = substruct('.', fields{i});
        val  = builtin('subsref', obj, sub1);
        sub2 = findfield(val, field);
        if numel(sub2) > 0
            S = [sub1 sub2];
            return
        end
    end
end
S = [];
