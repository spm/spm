function cfg_ui_restore(en)
%CFG_UI_RESTORE Restore state of properties
% CFG_UI_RESTORE(en) restores property values that were disabled by
% CFG_UI_DISABLE.
%
% See also CFG_UI_DISABLE.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_ui_restore.m 5768 2013-11-26 12:29:11Z volkmar $

rev = '$Rev: 5768 $'; %#ok<NASGU>
if ~isempty(en) && isstruct(en)
    sel = ishandle(en.c);
    set(en.c(sel), en.property, 'on');
end
