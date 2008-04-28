function varargout=subsref(this,subs)
% SUBSREF Subscripted reference
% An overloaded function...
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Stefan Kiebel
% $Id: subsref.m 1491 2008-04-28 16:46:35Z vladimir $

if isempty(subs)
    return;
end;

if this.Nsamples == 0
    error('Attempt to reference a field of an empty meeg object.');
end;

switch subs(1).type
    case '()'
        if numel(subs)~= 1, error('Expression too complicated');end;
        varargout = {double(subsref(this.data.y, subs))};
    case '{}'
    case '.'
        if isfield(this.other, subs(1).subs)
            field = getfield(this.other, subs(1).subs);
            if numel(subs)==1
                varargout = {field};
            else
                varargout ={subsref(field, subs(2:end))};
            end
        elseif ismethod(this, subs(1).subs)
            switch numel(subs)
                case 1
                    varargout = {feval(subs(1).subs, this)};
                case 2
                    varargout = {feval(subs(1).subs, this, subs(2).subs{:})};
                otherwise
                    error('Expression too complicated');
            end
        else
            error('Reference to non-existent or private meeg method or field.');
        end
    otherwise
        error('Unfamiliar referencing type');
end


