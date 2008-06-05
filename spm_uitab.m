function [handles] = spm_uitab(hparent,labels,callbacks,tag,active)
% function [handles] = spm_uitab(hfig,labels,callbacks)
% This functiuon creates tabs in the SPM graphics window.
% These tabs may be associated with different sets of axes and uicontrol,
% through the use of callback functions linked to the tabs.
% IN:
%   - hparent: the handle of the parent of the tabs (can be the SPM graphics
%   windows, or the handle of the uipanel of a former spm_uitab...)
%   - labels: a cell array of string containing the labels of the tabs
%   - callbacks: a cell array of strings which will be evaluated using the
%   'eval' function when clicking on a tab.
%   - tag: a string which is the tags associated with the tabs (useful for
%   finding them in a window...)
%   - active: the index of the active tab when creating the uitabs (default
%   = 1).
% OUT:
%   - handles: a structure of handles for the differents tab objects.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_uitab.m 1792 2008-06-05 15:43:54Z jean $

Ntabs = length(labels);

if ~exist('callbacks','var') || isempty(callbacks)
    for i=1:Ntabs
        callbacks{i} = [];
    end
end
if  ~exist('tag','var') || isempty(tag)
    tag = '';
end
if  ~exist('active','var') || isempty(active)
    active = 1;
end

if ~isequal(get(hparent,'type'),'figure')
    POS = get(hparent,'position');
    pos1 = [POS(1)+0.02,POS(2)+0.01,POS(3)-0.04,POS(4)-0.06];
else
    pos1 = [0.01 0.005 0.98 0.965];
end


handles.hp = uipanel('position',pos1,...
    'BorderType','beveledout',...
    'BackgroundColor',0.95*[1 1 1],...
    'parent',hparent,...
    'tag',tag);
set(handles.hp,'units','normalized');

xl = pos1(1);
yu = pos1(2) +pos1(4);
dx = 0.1;
ddx = 0.0025;
ddy = 0.005;
dy = 0.025;
for i =1:Ntabs
    handles.htab(i) = uicontrol('style','pushbutton',...
        'units','normalized','position',[xl+dx*(i-1) yu dx dy],...
        'SelectionHighlight','off',...
        'string',labels{i},...
        'BackgroundColor',0.95*[1 1 1],...
        'parent',hparent,...
        'tag',tag);
    set(handles.htab(i),'units','normalized')
    handles.hh(i) = uicontrol('style','text',...
        'units','normalized','position',...
        [xl+dx*(i-1)+ddx yu-ddy dx-2*ddx 2*ddy],...
        'BackgroundColor',0.95*[1 1 1],...
        'parent',hparent,...
        'tag',tag);
    set(handles.hh(i),'units','normalized')
end
set(handles.hh(active),'visible','on')
others = setdiff(1:Ntabs,active);
set(handles.htab(active),...
    'FontWeight','bold');
set(handles.hh(others),'visible','off');
set(handles.htab(others),...
    'ForegroundColor',0.25*[1 1 1]);
ud.handles = handles;
for i =1:Ntabs
    ud.ind = i;
    ud.callback = callbacks{i};
    set(handles.htab(i),'callback',@doChoose,'userdata',ud,...
        'BusyAction','cancel',...
        'Interruptible','off');
end



function doChoose(o1,o2)
ud = get(gco,'userdata');
% Do nothing if called tab is curret (active) tab
if ~strcmp(get(ud.handles.htab(ud.ind),'FontWeight'),'bold')
    set(ud.handles.hh(ud.ind),'visible','on');
    set(ud.handles.htab(ud.ind),...
        'ForegroundColor',0*[1 1 1],...
        'FontWeight','bold');
    others = setdiff(1:length(ud.handles.hh),ud.ind);
    set(ud.handles.hh(others),'visible','off');
    set(ud.handles.htab(others),...
        'ForegroundColor',0.25*[1 1 1],...
        'FontWeight','normal');
    if ~isempty(ud.callback)
        eval(ud.callback)
    end
    drawnow
end

