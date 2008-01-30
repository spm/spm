function obj = subsasgn(obj,subs,dat)
% Overloaded subsasgn function for meeg objects.
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: subsasgn.m 1125 2008-01-30 12:12:18Z vladimir $

if isempty(subs)
    return;
end;

if obj.Nsamples == 0
    error('Attempt to assign to a field of an empty meeg object.');
end;


if strcmp(subs(1).type, '.')
    if ismethod(obj, subs(1).subs)
        error('meeg method names cannot be used for custom fields');
    else
        if isempty(obj.other)
            obj.other = struct(subs(1).subs, dat);
        else
            obj.other = subsasgn(obj.other, subs, dat);
        end
    end
else
    error('Unsupported assignment type for meeg.');
end



