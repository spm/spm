function this = subsasgn(this,subs,dat)
% Overloaded subsasgn function for meeg objects.
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: subsasgn.m 4432 2011-08-15 12:43:44Z christophe $

if isempty(subs)
    return;
end;

if this.Nsamples == 0
    error('Attempt to assign to a field of an empty meeg object.');
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
    if numel(subs)~= 1, error('Expression too complicated');end;
    if this.montage.Mind>0
        error('Attempt to assign to a meeg object with online montage applied.');
    end;    
    this.data.y = subsasgn(this.data.y, subs, dat);
else
    error('Unsupported assignment type for meeg.');
end



