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
% $Id: spm_ov_browser.m 5585 2013-07-19 15:50:44Z guillaume $


cmd = lower(varargin{1});
switch cmd
    % Context menu and callbacks
    case 'context_menu'
        ret = uimenu(varargin{3}, ...
            'Label',    'Browse...', ...
            'Tag',      'orthviews_browser', ...
            'Callback', @browser_ui);
    case 'ui'
        if nargin == 1
            browser_ui;
        else
            browser(varargin{2:end});
        end
    case 'redraw'
        browser_redraw(varargin{2:end});
    otherwise
end


%==========================================================================
function browser_ui(hObj,event)

[f,sts] = spm_select([2 Inf],'image','Select images...');
if ~sts, return; end

if nargin
    Fgraph = ancestor(hObj,'figure');
    hC = current_handle;
else
    spm_check_registration(f(1,:));
    Fgraph = [];
    hC = 1;
end

hS = browser(f, Fgraph, hC);

hM = getappdata(hS,'hM');
set(hM,'Label','Browse','Callback','');
h = uimenu('Parent',hM,'Label','Display profile','Callback',@browser_profile);
setappdata(h,'hS',hS);
h = uimenu('Parent',hM,'Label','Quit','Callback',@browser_quit_button);
setappdata(h,'hS',hS);


%==========================================================================
function hS = browser(f, Fgraph, hC)

global st
f = cellstr(f);
if nargin < 2 || isempty(Fgraph)
    Fgraph = st.fig;
end
if nargin < 3 || isempty(hC)
    hC = 1;
end

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

setappdata(hS,'f', f);
setappdata(hS,'hT',hT);
setappdata(hS,'hC',hC);
setappdata(hS,'hB',hB);
setappdata(hS,'hP',hP);
setappdata(hS,'hM',findobj(st.fig,'Type','uimenu','Tag','orthviews_browser'));


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
global st
hS = getappdata(hObj,'hS');
hC = getappdata(hS,'hC');
delete(getappdata(hS,'hT'));
delete(getappdata(hS,'hB'));
delete(getappdata(hS,'hP'));
hM = getappdata(hS,'hM');
set(hM,'Label','Browse...','Callback',@browser_ui);
delete(get(hM,'Children'));
delete(hS);
try, st.vols{hC} = rmfield(st.vols{hC},'browser'); end % remove redraw callback


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
function browser_redraw(i,varargin) %i, TM0, TD, CM0, CD, SM0, SD
global st
feval(st.vols{i}.browser.fun,st.vols{i}.browser.h);


%==========================================================================
function browser_profile(hObj,event)
global st
hS  = getappdata(hObj,'hS');
hC  = getappdata(hS,'hC');
hV  = getappdata(hObj,'hV');
if isempty(hV)
    hV = spm_vol(char(getappdata(hS,'f')));
    setappdata(hObj,'hV',hV);
end

pos = spm_orthviews('pos',hC);
Y = spm_get_data(hV,[pos(1),pos(2),pos(3)]',false);

hAx = getappdata(hObj,'hAx');
if isempty(hAx) || ~ishandle(hAx)
    hF = figure;
    hAx = axes('Parent',hF);
    setappdata(hObj,'hAx',hAx);
end

plot(hAx,Y);
ylabel(hAx,sprintf('[%.2f %.2f %.2f]',pos));
st.vols{hC}.browser.fun = @browser_profile;
st.vols{hC}.browser.h = hObj;


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
