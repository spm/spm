function varargout=subsref(obj,subs)
% SUBSREF Subscripted reference
% An overloaded function...
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: subsref.m 1125 2008-01-30 12:12:18Z vladimir $

if isempty(subs)
    return;
end;

if obj.Nsamples == 0
    error('Attempt to reference a field of an empty meeg object.');
end;

switch subs(1).type
    case '()'
        if numel(subs)~= 1, error('Expression too complicated');end;
        varargout = {subsref(obj.data.y, subs)};
    case '{}'
    case '.'
        if ismethod(obj, subs(1).subs)
            switch numel(subs)
                case 1
                    varargout = {feval(subs(1).subs, obj)};
                case 2
                    varargout = {feval(subs(1).subs, obj, subs(2).subs{:})};
                otherwise
                    error('Expression too complicated');
            end
        elseif isfield(obj.other, subs(1).subs)
            field = getfield(obj.other, subs(1).subs);
            if numel(subs)==1
                varargout = {field};
            else
                varargout ={subsref(field, subs(2:end))};
            end
        else
            error('Reference to non-existent or private meeg method or field.');
        end
    otherwise
        error('Unfamiliar referencing type');
end


