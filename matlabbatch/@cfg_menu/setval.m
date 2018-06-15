function item = setval(item, val, dflag)

% function item = setval(item, val, dflag)
% set item.val{1} to item.values{val}. If val == {}, set item.val to {}
% If dflag is true, and item.cfg_item.def is not empty, set the default setting for
% this item instead by calling feval(item.cfg_item.def{:}, val). If val == {}, use
% the string '<UNDEFINED>' as in a harvested tree. If dflag is true, but
% no item.cfg_item.def defined, set item.val{1} instead.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: setval.m 7342 2018-06-15 12:44:44Z volkmar $

rev = '$Rev: 7342 $'; %#ok

if iscell(val) && isempty(val)
    if dflag
        if ~isempty(item.cfg_item.def)
            try
                feval(item.cfg_item.def, {'<UNDEFINED>'});
            catch
                cfg_message('matlabbatch:setval:defaults', ...
                            '%s: unable to set default value.', ...
                            subsasgn_checkstr(item, substruct('.','val')));
            end
        else
            item = subsasgn(item, substruct('.','val'), {});
        end
    else
        item = subsasgn(item, substruct('.','val'), {});
    end
else
    if ~isa(val, 'cfg_dep')
        if isnumeric(val) && val >= 1 && val <= numel(item.values)
            val = item.values{val};
        else
            cfg_message('matlabbatch:setval', ...
                '%: unable to set value. Index must be an integer between 1 and %d.', ...
                subsasgn_checkstr(item, substruct('.','val')), numel(item.values));
            return
        end
    end
    if dflag
        [sts, val1] = subsasgn_check(item, substruct('.','val'), {val});
        if sts
            if ~isempty(item.cfg_item.def)
                try
                    feval(item.cfg_item.def, val1);
                catch
                    cfg_message('matlabbatch:setval:defaults', ...
                                '%s: unable to set default value.', ...
                                subsasgn_checkstr(item, substruct('.','val')));
                end
            else
                item = subsasgn(item, substruct('.','val'), {val});
            end
        end
    else
        item = subsasgn(item, substruct('.','val'), {val});
    end
end
