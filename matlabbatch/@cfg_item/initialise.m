function item = initialise(item, val, dflag)

% function item = initialise(item, val, dflag)
% This is a generic initialisation function to insert values into the val
% field of a cfg_item. It is usable for all leaf cfg_item classes
% (cfg_entry, cfg_files, cfg_menu). Assignment checks are done through
% subsasgn.
% dflag is ignored for leaf entries - the passed value is assigned
% unconditionally.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: initialise.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

subs = substruct('.', 'val', '{}', {1});
if ischar(val) && strcmp(val, '<UNDEFINED>') % val may be <UNDEFINED>
    item.val = {};
else
    item = subsasgn(item, subs, val);
end;