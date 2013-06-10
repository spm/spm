function ret = spm_ov_browser(varargin)
% Browser tool - plugin for spm_orthviews
%
% This routine is a plugin to spm_orthviews. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the MATLAB prompt.
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_ov_browser.m 5534 2013-06-10 13:29:09Z guillaume $


cmd = lower(varargin{1});
switch cmd
    % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, ...
            'Label', 'Browse...', ...
            'Callback', @browser);
        ret = item0;
    otherwise
end


%==========================================================================
function browser(hObj,event)

[f,sts] = spm_select([2 Inf],'image','Select images...');
if ~sts, return; end
f = cellstr(f);

Fgraph = ancestor(hObj,'figure');

hS = uicontrol('Parent', Fgraph,...
    'Style',             'slider',...
    'Units',             'normalized',...
    'Position',          [0.1 0.025 0.8 0.02],...
    'Min',               1,...
    'Max',               numel(f),...
    'Value',             1,...
    'SliderStep',        [1 1]/(numel(f)-1));
try
    hListener = handle.listener(hS,'ActionEvent',@browser_slider);
    setappdata(hS,'myListener',hListener);
catch
    set(hS,'Callback',   @browser_slider);
end
    
hB = uicontrol('Parent', Fgraph,...
    'Style',             'togglebutton',...
    'Units',             'normalized',...
    'Position',          [0.92 0.025 0.03 0.02],...
    'String',            '>',...
    'Callback',          @browser_play_button);
setappdata(hB,'hS',hS);

hP = uicontrol('Parent', Fgraph,...
    'Style',             'pushbutton',...
    'Units',             'normalized',...
    'Position',          [0.96 0.025 0.03 0.02],...
    'String',            'Q',...
    'TooltipString',     'Quit',...
    'Callback',          @browser_quit_button);
setappdata(hP,'hS',hS);

hT = uicontrol('Parent',   Fgraph,...
    'Style',               'text',...
    'Units',               'normalized',...
    'Position',            [0.1 0.0025 0.8 0.02],...
    'HorizontalAlignment', 'center',...
    'BackgroundColor',     [1 1 1],...
    'String',              f{1});

hC = current_handle;

setappdata(hS,'f', f);
setappdata(hS,'hT',hT);
setappdata(hS,'hC',hC);
setappdata(hS,'hB',hB);
setappdata(hS,'hP',hP);


%==========================================================================
function browser_play_button(hObj,event)
hS = getappdata(hObj,'hS');
f  = getappdata(hS,'f');
j  = round(get(hS,'Value'));
if j == numel(f), j = 1; end % if at end already, play from start
tp = 1 / numel(f); % make the complete sequence take at least 1 second
for i=j:numel(f)
    if ~get(hObj,'Value'), return; end
    t = tic;
    set(hS,'Value',i);
    browser_slider(hS);
    pause(tp - toc(t))
end
set(hObj,'Value',0);


%==========================================================================
function browser_quit_button(hObj,event)
hS = getappdata(hObj,'hS');
delete(getappdata(hS,'hT'));
delete(getappdata(hS,'hB'));
delete(getappdata(hS,'hP'));
delete(hS);

%==========================================================================
function browser_slider(hObj,event)
global st
i  = round(get(hObj,'Value'));
f  = getappdata(hObj,'f');
hT = getappdata(hObj,'hT');
hC = getappdata(hObj,'hC');

V  = spm_vol(f{i});
fn = fieldnames(V);
for k=1:numel(fn)
    st.vols{hC}.(fn{k}) = V.(fn{k});
end
hM = findobj(st.vols{hC}.ax{1}.cm,'UserData','filename');
spm_orthviews('context_menu','image_info',hM,hC);
pos = spm_orthviews('pos',hC);
set(findobj(st.vols{hC}.ax{1}.cm,'UserData','v_value'),...
    'Label',sprintf('Y = %g',spm_sample_vol(st.vols{hC},pos(1),pos(2),pos(3),st.hld)));
spm_orthviews('Redraw');
set(hT,'String',f{i});


%==========================================================================
function h = current_handle
global st
try
    hs = [];
    for i = spm_orthviews('valid_handles')
        hs = [hs st.vols{i}.ax{1}.cm];
    end
    hc = get(gca,'UIContextMenu');
    h  = find(hs==hc);
catch
    h  = 1;
end
