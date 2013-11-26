function en = cfg_ui_disable(hObject, property)
%CFG_UI_DISABLE Disable properties
% en = CFG_UI_DISABLE(hObject, property) disables property in all children
% of hObject, returning their handles in en.c and previous state in cell
% list en.en. CFG_UI_RESTORE(en) can be used to restore the property to
% their original setting.
% Property must be a property that has the values 'on' (enabled) or 'off'
% (disabled).
%
% See also CFG_UI_RESTORE.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_ui_disable.m 5768 2013-11-26 12:29:11Z volkmar $

rev = '$Rev: 5768 $';  %#ok<NASGU>
en.c        = findall(hObject, property, 'on');
en.property = property;
set(en.c, en.property, 'off');

