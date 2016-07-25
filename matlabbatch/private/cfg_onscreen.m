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
% $Id: cfg_onscreen.m 6840 2016-07-25 12:21:25Z guillaume $

rev = '$Rev: 6840 $'; 

% save figure units - use pixels here
units = get(fg,'Units');
set(fg,'Units','pixels');
Rect = get(fg,'Position');
S0   = get(0,'MonitorPositions');
if size(S0,1) > 1 % Multiple Monitors
    %-Use Monitor containing the Pointer
    pl = get(0,'PointerLocation');
    w  = find(pl(1)>=S0(:,1) & pl(1)<S0(:,1)+S0(:,3)-1 &...
            pl(2)>=S0(:,2) & pl(2)<S0(:,2)+S0(:,4));
    if numel(w)~=1, w = 1; end
    S0 = S0(w,:);
end
if ~(Rect(1)>=S0(1) && Rect(1)+Rect(3)<=S0(3) && Rect(2)>=S0(2) && Rect(2)+Rect(4)<=S0(4))
    Rect(1) = S0(1) + (S0(3) - Rect(3))/2;
    Rect(2) = S0(2) + (S0(4) - Rect(4))/2;
end
if nargout == 0
    set(fg, 'Position',Rect);
    figure(fg);
else
    varargout{1} = Rect;
end
set(fg,'Units',units);
