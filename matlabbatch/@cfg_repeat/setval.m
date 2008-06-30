function item = setval(item, val, dflag)

% function item = setval(item, val, dflag)
% Add, replicate or delete an entry in item.val. The semantics is based on
% the contents of the 2nd argument:
% If val == {}, set item.val to {}.
% If val(1) > 0, set item.val{min(val(2), numel(item.val)+1)} to
% item.values{val(1)}. 
% If val(1) =< 0, replicate item.val{min(val(2), numel(item.val))} by
% appending it to the list in item.val. Note that no provision is being
% made to clear fields in cfg_exbranches that deal with dependencies and
% the overall job structure. This has to be done in the job management
% utility.
% If val(1) is not finite, then the entry val(2) is deleted from item.val.
% dflag is ignored for cfg_repeat items.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: setval.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

if iscell(val) && isempty(val)
    item = subsasgn(item, substruct('.','val'), {});
else
    if val(1) > 0 && isfinite(val(1))
        val1 = item.values{val(1)};
        substype = '{}';
        % Add after last element, if out of range of val list
        subsind  = min(val(2),numel(subsref(item, substruct('.', 'val')))+1);
    elseif val(1) <= 0 && isfinite(val(1))
        % Get element to replicate
        replind = min(val(2),numel(subsref(item, substruct('.', 'val'))));
        if replind > 0 % there may be nothing to replicate
            val1 = item.cfg_item.val{replind};
            % Add after last element
            substype = '{}';
            subsind = numel(subsref(item, substruct('.', 'val')))+1;
        else
            cfg_message('matlabbatch:setval:repl','Nothing to replicate.');
            return;
        end;
    else
        val1 = [];
        substype = '()';
        % Delete last element, if out of range of val list
        subsind  = min(val(2),numel(subsref(item, substruct('.', 'val'))));
    end
    item = subsasgn(item, substruct('.','val', substype, {subsind}), val1);
end;
