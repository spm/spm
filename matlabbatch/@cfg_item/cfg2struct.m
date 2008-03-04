function sitem = cfg2struct(item)

% function sitem = cfg2struct(item)
% Return a struct containing all fields of item plus a field type. This is
% the generic method.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg2struct.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

% Remove class property, save it as field 'type'
sitem = struct(item);
sitem.type = class(item);

% Treat val{:} fields, if they are cfg_item objects
% Assume that all of them are cfg_items, if the first one is
if numel(item.val) > 0 && isa(item.val{1}, 'cfg_item')
    for k = 1:numel(item.val)
        sitem.val{k} = cfg2struct(item.val{k});
    end;
end;
