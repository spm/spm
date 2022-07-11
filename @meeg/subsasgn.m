function this = subsasgn(this,subs,dat)
% Overloaded subsasgn function for meeg objects.
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if isempty(subs)
    return;
end;

if strcmp(subs(1).type, '.')
    if ismethod(this, subs(1).subs)
        error('meeg method names cannot be used for custom fields');
    else
        if isempty(this.other)
            this.other = struct(subs(1).subs, []);
            this.other = builtin('subsasgn', this.other, subs, dat);
        else
            this.other = builtin('subsasgn', this.other, subs, dat);
        end
    end
elseif strcmp(subs(1).type, '()')
    if ~islinked(this), error('The object is not linked to data file'); end
    if numel(subs)~= 1,  error('Expression too complicated');           end
    if this.montage.Mind>0
        error('Attempt to assign to a meeg object with online montage applied.');
    end
    this.data = subsasgn(this.data, subs, dat);
else
    error('Unsupported assignment type for meeg.');
end
