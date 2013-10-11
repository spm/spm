function [item, inputs] = fillvals(item, inputs, infcn)

% function [item, inputs] = fillvals(item, inputs, infcn)
% If ~all_set_item, try to set item.val to the items listed in inputs{1}.
% inputs{1} should be a cell array of indices into item.values. For
% cfg_choice items, this list should only contain one item.
% Validity checks are performed through setval. If inputs{1} is not
% suitable for this item, it is discarded. If infcn is a function handle,
% [val sts] = infcn(item) 
% will be called to obtain a value for this item. This call will be
% repeated until either val can be assigned to item or sts is true.
%
% This function is identical for all cfg_intree classes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: fillvals.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

% Set item itself
if ~all_set_item(item)
    if ~isempty(inputs)
        for k = 1:numel(inputs{1})
            item = setval(item, [inputs{1}{k} Inf], false);
        end;
        inputs = inputs(2:end);
    end;
    if ~all_set_item(item) && ~isempty(infcn) && subsasgn_check_funhandle(infcn)
        sts = false;
        while ~sts && ~all_set_item(item)
            [val, sts] = feval(infcn, item);
            if sts
                for k = 1:numel(val)
                    item = setval(item, [val{k} Inf], false);
                end;
            end;
        end;
    end;
end;
% Set children
citems = subsref(item, substruct('.','val'));
for k = 1:numel(citems)
    if ~all_set(citems{k})
        [citems{k}, inputs] = fillvals(citems{k}, inputs, infcn);
    end;
end;
item = subsasgn(item, substruct('.','val'), citems);
