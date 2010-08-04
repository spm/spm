function varargout = cfg_onscreen(fg)
% Move figure on the screen containing the mouse
%    cfg_onscreen(fg) - move figure fg on the screen containing the mouse
%    pos = cfg_onscreen(fg) - compute position of figure, do not move it
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_onscreen.m 4033 2010-08-04 15:53:35Z volkmar $

rev = '$Rev: 4033 $'; 

% save figure units - use pixels here
units = get(fg,'Units');
set(fg,'Units','pixels');
Rect = get(fg,'Position');
S0   = get(0,'MonitorPosition');
if size(S0,1) > 1 % Multiple Monitors
    %-Use Monitor containing the Pointer
    pl = get(0,'PointerLocation');
    w  = find(pl(1)>=S0(:,1) & pl(1)<S0(:,1)+S0(:,3)-1 &...
            pl(2)>=S0(:,2) & pl(2)<S0(:,2)+S0(:,4));
    if numel(w)~=1, w = 1; end
    S0 = S0(w,:);
end
Rect(1) = S0(1) + (S0(3) - Rect(3))/2;
Rect(2) = S0(2) + (S0(4) - Rect(4))/2;
if nargout == 0
    set(fg, 'Position',Rect);
else
    varargout{1} = Rect;
end
set(fg,'Units',units);
