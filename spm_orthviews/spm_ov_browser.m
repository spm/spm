function ret = spm_ov_browser(varargin)
% Browser tool - plugin for spm_orthviews
%
% This routine is a plugin to spm_orthviews. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the MATLAB prompt.
%__________________________________________________________________________
% Copyright (C) 2013-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_ov_browser.m 6750 2016-03-21 12:36:08Z guillaume $


if ~nargin, varargin = {'ui'}; end
cmd = lower(varargin{1});
switch cmd
    % Context menu and callbacks
    case 'context_menu'
        ret = uimenu(varargin{3}, ...
            'Label',    'Browse...', ...
            'Tag',      'orthviews_browser', ...
            'Callback', @browser_ui);
    case 'ui'
        if nargin <= 1
            browser_ui;
        else
            browser_ui([],[],varargin{2:end});
        end
    case 'redraw'
        browser_redraw(varargin{2:end});
    otherwise
end


%==========================================================================
function browser_ui(hObj,event,varargin)

global st
if nargin < 3
    if nargin && ~isempty(hObj)
        hC = current_handle;
        sel = {[st.vols{hC}.fname ',' num2str(st.vols{hC}.n(1))]};
    else
        sel = {};
    end
    [f,sts] = spm_select([2 Inf],'image','Select images...',sel);
    if ~sts, return; end
    f = cellstr(f);
else
    f = varargin;
    if numel(f) == 1
        if ischar(f{1})
            f = cellstr(f{1});
        elseif iscellstr(f{1})
            f = f{1};
        end
    end
end

if nargin && ~isempty(hObj)
    Fgraph = ancestor(hObj,'figure');
    hC = current_handle;
else
    spm_check_registration(f{1});
    Fgraph = [];
    hC = 1;
end

hS = browser(f, Fgraph, hC);

hM = getappdata(hS,'hM');
for i=1:numel(hM)
    set(hM,'Label','Browse','Callback','');
    h = uimenu('Parent',hM(i),'Label','Display profile','Callback',@browser_profile);
    setappdata(h,'hS',hS);
    h = uimenu('Parent',hM(i),'Label','Save movie...','Callback',@browser_movie);
    setappdata(h,'hS',hS);
    h = uimenu('Parent',hM(i),'Label','Loop','Checked','off','Callback',@browser_loop);
    setappdata(h,'hS',hS);
    h = uimenu('Parent',hM(i),'Label','Quit','Callback',@browser_quit_button);
    setappdata(h,'hS',hS);
end


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
    if spm_check_version('matlab','8.4') >= 0
        hListener = addlistener(hS,'ContinuousValueChange',@browser_slider);
    else
        hListener = handle.listener(hS,'ActionEvent',@browser_slider);
    end
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
n  = numel(f);
if j == n, j = 1; end % if at end already, play from start
tp = 1 / n; % make the complete sequence take at least 1 second
hM = getappdata(hS,'hM');
hL = findobj(hM,'Label','Loop');
if strcmp(get(hL,'Checked'),'on'), l = intmax; else l = n; end
for i=j:l
    try, if ~get(hObj,'Value'), return; end; catch, return; end
    t = tic;
    k = mod(i,n);
    if k == 0, k = n; end
    set(hS,'Value',k);
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
for i=1:numel(hM)
    set(hM(i),'Label','Browse...','Callback',@browser_ui);
    delete(get(hM(i),'Children'));
end
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
try, set(hT,'String',f{i}); end


%==========================================================================
function browser_redraw(i,varargin) %i, TM0, TD, CM0, CD, SM0, SD
global st
feval(st.vols{i}.browser.fun,st.vols{i}.browser.h);


%==========================================================================
function browser_loop(hObj,event)
if strcmpi(get(hObj,'Checked'),'off')
    set(hObj,'Checked','on');
else
    set(hObj,'Checked','off');
end


%==========================================================================
function browser_profile(hObj,event)
global st
if nargin > 1
    if strcmpi(get(hObj,'Checked'),'off')
        set(hObj,'Checked','on');
    else
        set(hObj,'Checked','off');
        hS = getappdata(hObj,'hS');
        hC = getappdata(hS,'hC');
        try, st.vols{hC} = rmfield(st.vols{hC},'browser'); end % remove redraw callback
        return;
    end
end
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
    %hM = uimenu('Parent',hF,'Label','&Options');
    %hM = uimenu('Parent',hM,'Label','Select SPM.mat...','Callback',@browser_profile_callback);
    %setappdata(hM,'hS',hS);
    %setappdata(hM,'hAx',hAx);
end

l = plot(hAx,Y);
set(l,'Tag','profile');
hold(hAx,'on')
i = round(get(hS,'value'));
plot(hAx,i,Y(i),'r*');
hold(hAx,'off');
ylabel(hAx,sprintf('[%.2f %.2f %.2f]',pos));
st.vols{hC}.browser.fun = @browser_profile;
st.vols{hC}.browser.h = hObj;


%==========================================================================
function browser_movie(hObj,event)
global st
hS = getappdata(hObj,'hS');
hC = getappdata(hS,'hC');
hS = getappdata(hObj,'hS');
f  = getappdata(hS,'f');

[filename, pathname] = uiputfile(...
    {'*.gif' 'GIF files (*.gif)';...
    '*.png' 'PNG files (*.png)'; ...
    '*.avi' 'AVI files (*.avi)';...
    }, 'Save as');
if isequal(filename,0) || isequal(pathname,0), return; end

file = fullfile(pathname,filename);

p1 = get(st.vols{hC}.ax{1}.ax,'Position');
p2 = get(st.vols{hC}.ax{3}.ax,'Position');
a  = [p1(1) p1(2)  p2(1)+p2(3)-p1(1) p2(2)+p2(4)-p1(2)] + 0.005*[-1 -1 2 2];
a  = max(min(a,1),0);

if strcmp(spm_file(file,'ext'),'avi')
    outputtype = 1;
    writerObj  = VideoWriter(file);
    open(writerObj);
elseif strcmp(spm_file(file,'ext'),'gif')
    outputtype = 2;
    delay = spm_input('Delay','+0','r',0.05,1,[0,Inf]);
    loop = spm_input('Loop','+1','r',0,1,[0,Inf]);
elseif strcmp(spm_file(file,'ext'),'png')
    outputtype = 3;
else
    error('Unknown output file type.');
end

for i=1:numel(f)
    set(hS,'Value',i);
    browser_slider(hS);
    
    X  = frame2im(getframe(st.fig));
    sz = size(X);
    sz = [sz(1) sz(1) sz(2) sz(2)];
    sz = ([1-a(2)-a(4),1-a(2),a(1),a(1)+a(3)] .* (sz-1)) + 1;
    sz = round(sz);
    X  = X(sz(1):sz(2),sz(3):sz(4),:);
    
    if outputtype == 1
        writeVideo(writerObj,X);
    elseif outputtype ==2
        [A,map] = rgb2ind(X,256);
        if i == 1
            imwrite(A,map,file,'gif','DelayTime',delay,'LoopCount',loop);
        else
            imwrite(A,map,file,'gif','DelayTime',delay,'WriteMode','append');
        end
    elseif outputtype ==3
        imwrite(X,spm_file(file,'suffix',sprintf('_%04d',i)),'png');
    end
end

if outputtype == 1, close(writerObj); end

fprintf('Movie saved in %s\n',file);                                    %-#


%==========================================================================
function browser_profile_callback(hObj,event)
[SPM,sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat...');
load(SPM);
if ~sts, return; end
hS = getappdata(hObj,'hS');
hC = getappdata(hS,'hC');
f  = getappdata(hS,'f');
hAx = getappdata(hObj,'hAx');

Y = get(findobj(hAx,'Tag','profile'),'YData')';
Y = spm_filter(SPM.xX.K,Y);
Y = SPM.xX.X*(pinv(spm_filter(SPM.xX.K,SPM.xX.X))*Y);
hold(hAx,'on');
plot(Y,'g');

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
