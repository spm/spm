function en = cfg_ui_disable(hObject, property)
%CFG_UI_DISABLE Disable properties
% en = CFG_UI_DISABLE(hObject, property) disables property in all children
% of hObject, returning their handles in en.c and previous state in cell
% list en.en. CFG_UI_RESTORE(en) can be used to restore the property to
% their original setting.
%
% See also CFG_UI_RESTORE.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_ui_disable.m 5679 2013-10-11 14:58:14Z volkmar $

rev = '$Rev: 5679 $';  %#ok<NASGU>
c   = findall(hObject);
sel = isprop(c,property);
en.c        = c(sel);
en.property = property;
en.en       = get(en.c, en.property);
set(en.c, en.property,'off');

