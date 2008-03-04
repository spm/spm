function item = initialise(item, val, dflag)

% function item = initialise(item, val, dflag)
% Disassemble val struct and pass its fields on to child nodes. If dflag
% is set, then val fields are passed on to the item.values{:} cfg_item(s),
% whose tags match a field name of val. In this case, multiple
% item.values{:} cfg_items can be initialised, if val has more than one
% field.
% If dflags is not set, then item.val is set to the initialised copy of the
% matching values item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: initialise.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

vtags = fieldnames(val);

if dflag % set defaults
    for k = 1:numel(item.values)
        % find field in val that corresponds to one of the branch vals
        vi = strcmp(gettag(item.values{k}), vtags);
        if any(vi) % field names are unique, so there will be at most one match
            item.values{k} = initialise(item.values{k}, ...
                val.(vtags{vi}), dflag);
        end;
    end;
else
    % select matching values struct, initialise and assign it to val
    % field
    for k = 1:numel(item.values)
        if strcmp(gettag(item.values{k}), vtags{1})
            item.cfg_item.val{1} = initialise(item.values{k}, ...
                val.(vtags{1}), dflag);
            break;
        end;
    end;
end;
