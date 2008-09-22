function ext = cfg_maxextent(obj, str)
% CFG_MAXEXTENT Returns the maximum extent of cellstr STR 
% Returns the maximum extent of obj OBJ when the cellstr STR will be
% rendered in it. MATLAB is not able to work this out correctly on its own
% for multiline strings. Therefore each line will be tried separately and
% its extent will be returned. To avoid 'flicker' appearance, OBJ should be
% invisible. The extent does not include the width of a scrollbar.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_maxextent.m 2131 2008-09-22 06:04:53Z volkmar $

rev = '$Rev: 2131 $'; %#ok
ext = zeros(size(str));
for k = 1:numel(str)
    set(obj,'String',str{k});
    next = get(obj, 'Extent');
    ext(k) = next(3);
end;
