function [ud] = spm_DisplayTimeSeries(y,options)
% Build a GUI for 'smart' time series display
% FORMAT [ud] = spm_DisplayTimeSeries(y,options)
% IN:
%   - y: the txn data, where t is the number of time sample, and p the
%        number of 'channels'
%   - options: a structure (default is empty), which allows to adapt this
%              function to specific needs. Optional fields are:
%       .hp: the handle of the parent figure/object. This is used to
%       include the time series display in a panel/figure. By default, a
%       new figure will be created.
%       .Fsample: the sample rate of the data (in Hz)
%       .events: a nex1 structure vector containing the time indices of the
%       events and their type (if any). Default is empty. Basic structure
%       contains fields .time and .type (see bellow).
%       .M: a pxn matrix premultiplied to the data when plotted (default is
%       1).
%       .bad a px1 binary vector containing the good/bad status of the
%       channels. Default is zeros(p,1).
%       .transpose: a binary variable that transposes the data (useful for
%       file_array display). Using options.transpose = 1 is similar to do
%       something similar to plot(y'). Default is 0.
%       .minY: the min value of the plotted data (used to define the main
%       axes limit). Default is calculated according to the offset.
%       .maxY: the max value of the plotted data (used to define the main
%       axes limit). Default is calculated according to the offset.
%       .minSizeWindow: minimum size of the plotted window (in number of
%       time samples). {min([200,0.5*size(y,1)]}
%       .maxSizeWindow: maximum size of the plotted window (in number of
%       time samples). {min([5000,size(y,1)])}
%       .ds: an integer giving the number of displayed time samples when
%       dynamically moving the display time window. Default is 1e4. If you
%       set it to Inf, no downsampling is applied.
%       .callback: a string or function handle which is evaluated after
%       each release of the mouse button (when moving the patch or clicking
%       on the slider). Default is empty.
%       .tag: a string used to tag both axes
%       .pos1: a 4x1 vector containing the position of the main display
%       axes {[0.13 0.3 0.775 0.655]}
%       .pos2: a 4x1 vector containing the position of the global power
%       display axes {[0.13 0.05 0.775 0.15]}
%       .pos3: a 4x1 vector containing the position of the temporal slider
%       {[0.13 0.01 0.775 0.02]}
%       .itw: a vector containing the indices of data time samples
%       initially displayed in the main axes {1:minSizeWindow}
%       .ytick: the 'ytick' property of the main axes
%       .yticklabel: the 'yticklabel' property of the main axes
%       .offset: a px1 vector containing the vertical offset that has to be
%       added to each of the plotted time series
%       !! .ytick, .yticklabel and .offset can be used to display labelled
%       time series one above each other !!
% OUT:
%   - ud: a structure containing all relevant informations about the
%   graphical objects created for the GUI. This is useful for manipulating
%   the figure later on (see below).
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_DisplayTimeSeries.m 6278 2014-12-04 13:20:20Z guillaume $


if ~exist('options','var')
    options = [];
end

if ~isa(y,'meeg')
    y(isnan(y)) = 0;
end

% Get optional parameters (if any) and set up defaults
%==========================================================================

if ~isempty(options) && isfield(options,'hp')
    hp = options.hp;
else
    hp = figure;
end

if ~isempty(options) && isfield(options,'Fsample')
    Fsample = options.Fsample;
else
    Fsample = 1;
end

if ~isempty(options) && isfield(options,'timeOnset')
    timeOnset = options.timeOnset;
else
    timeOnset = 0;
end

if ~isempty(options) && isfield(options,'events')
    events = options.events;
else
    events = [];
end

if ~isempty(options) && isfield(options,'transpose')
    transpose = options.transpose;
else
    transpose = 0;
end

if ~isempty(options) && isfield(options,'M')
    M = options.M;
    nc = size(M,1);
    if ~transpose
        nt = size(y,1);
    else
        nt = size(y,2);
    end
else
    M = 1;
    if ~transpose
        [nt,nc] = size(y);
    else
        [nc,nt] = size(y);
    end
end

if ~isempty(options) && isfield(options,'bad')
    bad = options.bad;
else
    bad = zeros(nc,1);
end

if ~isempty(options) && isfield(options,'minSizeWindow')
    minSizeWindow = options.minSizeWindow;
else
    minSizeWindow = 200;
end
minSizeWindow = min([minSizeWindow,.5*nt]);

if ~isempty(options) && isfield(options,'maxSizeWindow')
    maxSizeWindow = options.maxSizeWindow;
else
    maxSizeWindow = 5000;
end
maxSizeWindow = min([maxSizeWindow,nt]);

if ~isempty(options) && isfield(options,'ds')
    ds = options.ds;
else
    ds = 1e4;
end

if ~isempty(options) && isfield(options,'callback')
    callback = options.callback;
else
    callback = [];
end

if ~isempty(options) && isfield(options,'tag')
    tag = options.tag;
else
    tag = '';
end

if ~isempty(options) && isfield(options,'pos1')
    pos1 = options.pos1;
else
    pos1 = [0.13 0.32 0.775 0.655];
end

if ~isempty(options) && isfield(options,'pos2')
    pos2 = options.pos2;
else
    pos2 = [0.13 0.12 0.775 0.15];
end

if ~isempty(options) && isfield(options,'pos3')
    pos3 = options.pos3;
else
    pos3 = [0.10 0.02 0.84 0.05];
end

if ~isempty(options) && isfield(options,'itw')
    itw = options.itw;
else
    itw = 1:minSizeWindow;
end

if ~isempty(options) && isfield(options,'ytick')
    ytick = options.ytick;
else
    ytick = [];
end

if ~isempty(options) && isfield(options,'yticklabel')
    yticklabel = options.yticklabel;
else
    yticklabel = [];
end

if ~isempty(options) && isfield(options,'offset')
    offset = options.offset(:);
else
    offset = zeros(nc,1);
end

if ~isempty(options) && isfield(options,'minY')
    minY = options.minY;
else
    yo = y + repmat(offset',nt,1);
    minY = min(yo(:));
end

if ~isempty(options) && isfield(options,'maxY')
    maxY = options.maxY;
else
    try
        maxY = max(yo(:));
    catch
        yo = y + repmat(offset',nt,1);
        maxY = max(yo(:));
    end

end



% Initialize display
%==========================================================================

% Get basic info about data
ud.y    = y;
ud.v.bad = bad;
ud.v.transpose = transpose;
ud.v.nc = nc;
ud.v.nt = nt;
ud.v.ds = ds;
ud.v.M = M;
ud.v.mi = minY;
ud.v.ma = maxY;
ud.v.minSizeWindow = minSizeWindow;
ud.v.maxSizeWindow = maxSizeWindow;
ud.v.offset = offset;
ud.v.et = events;
ud.v.handles.hp = hp;
ud.callback = callback;

% Get downsampled global power
decim = max([1,round(nt./1000)]);
ud.v.ind = 1:decim:nt;
if ~ud.v.transpose
    My = ud.v.M*y(ud.v.ind,:)';
    ud.v.y2 = sum(My.^2,1);
else
    My = ud.v.M*y(:,ud.v.ind);
    ud.v.y2 = sum(My.^2,1);
end


mi = min(ud.v.y2);
ma = max(ud.v.y2);

if mi == 0 && ma == 0
    mi = -eps;
    ma = eps;
else
    mi = mi - mi.*1e-3;
    ma = ma + ma.*1e-3;
end

if spm_check_version('matlab','8.4') >= 0
    dispmode = {'SortMethod','childorder'};
else
    dispmode = {'drawmode','fast'};
end

% Create axes
ud.v.handles.axes = axes('parent',hp,...
    'units','normalized',...
    'position',pos1,...
    'xtick',[],'xticklabel',[],...
    'ytick',ytick,'yticklabel',yticklabel,...
    'tag',tag,...
    'nextplot','add',...
    dispmode{:});
ud.v.handles.gpa  = axes('parent',hp,...
    'units','normalized',...
    'position',pos2,...
    'tag',tag,'nextplot','add','ytick',[],...
    'box','off',...
    'color','none',...
    'ygrid','off',...
    dispmode{:});

% Initialize time series
col = colormap(lines);
col = col(1:7,:);
if ~ud.v.transpose
    My = ud.v.M*y(itw,:)';
else
    My = ud.v.M*y(:,itw);
end
for i=1:ud.v.nc
    ii = mod(i+7,7)+1;
    ud.v.handles.hp(i) = plot(ud.v.handles.axes,...
        itw,My(i,:)+offset(i),...
        'color',col(ii,:));
    if ud.v.bad(i)
        set(ud.v.handles.hp(i),'linestyle',':')
    end
end
ud.v.handles.gpp = plot(ud.v.handles.gpa,...
    ud.v.ind,ud.v.y2,...
    'color',0.5*[1 1 1]);


% Add events if any
if ~isempty(ud.v.et)
    for i=1:length(ud.v.et)
        ud.v.et(i).col = mod(ud.v.et(i).type+7,7)+1;
        ud.v.et(i).hp = plot(ud.v.handles.axes,...
            ud.v.et(i).time.*[1 1],...
            [ud.v.mi,ud.v.ma],...
            'color',col(ud.v.et(i).col,:),...
            'userdata',i,...
            'ButtonDownFcn','set(gco,''selected'',''on'')',...
            'Clipping','on');
    end
end

% Add display scrolling patch and borders
ud.v.handles.pa = patch('parent',ud.v.handles.gpa,...
    'xdata',[itw(1) itw(end) itw(end) itw(1)],...
    'ydata',[mi mi ma ma],...
    'edgecolor','none',...
    'facecolor',[.5 .5 .5],...
    'facealpha',0.5,...
    'ButtonDownFcn',@doPatch,...
    'interruptible','off');
ud.v.handles.lb = plot(ud.v.handles.gpa,...
    [itw(1) itw(1)],[mi ma],...
    'k',...
    'buttondownfcn',@doLb,...
    'interruptible','off');
ud.v.handles.rb = plot(ud.v.handles.gpa,...
    [itw(end) itw(end)],[mi ma],...
    'k',...
    'buttondownfcn',@doRb,...
    'interruptible','off');

% Adapt axes properties to display
tgrid = (0:(nt-1))./Fsample + timeOnset;
set(ud.v.handles.gpa,...
    'ylim',[mi ma],...
    'xlim',[1 nt]);
xtick = floor(get(ud.v.handles.gpa,'xtick'));
xtick(xtick==0) = 1;
a = cell(length(xtick),1);
for i=1:length(xtick)
    a{i} = num2str(tgrid(xtick(i)),2);
end
set(ud.v.handles.gpa,...
    'xtick',xtick,...
    'xticklabel',a);
set(ud.v.handles.axes,...
    'xlim',[itw(1) itw(end)],...
    'ylim',[ud.v.mi ud.v.ma]);

% Add temporal slider
ud.v.handles.hslider = uicontrol('parent',hp,...
    'style','slider',...
    'units','normalized',...
    'Position',pos3,...
    'min',max([1,minSizeWindow/2-1]),...
    'max',min([ud.v.nt,ud.v.nt-minSizeWindow/2+1]),...
    'value',mean([itw(1),itw(end)]),...
    'sliderstep',.1*[minSizeWindow/(ud.v.nt-1) 4*minSizeWindow/(ud.v.nt-1)],...
    'callback',@doScroll,...
    'userdata',ud.v.handles.gpa,...
    'BusyAction','cancel',...
    'Interruptible','on',...
    'tooltipstring','Scroll data',...
    'tag',tag);

% Store required info in global power axes
set(ud.v.handles.gpa,'userdata',ud)
axes(ud.v.handles.gpa)





% Subfunctions
%==========================================================================

function doScroll(src,evt)
gpa = get(src,'userdata');
xm = get(src,'value');
ud = get(gpa,'userdata');
sw = diff(get(ud.v.handles.axes,'xlim'));
xl = xm + [-sw./2,+sw./2];
if xl(1) >= 1 && xl(2) <= ud.v.nt
    xl = round(xl);
    if ~ud.v.transpose
        My = ud.v.M*ud.y(xl(1):xl(2),:)';
    else
        My = ud.v.M*ud.y(:,xl(1):xl(2));
    end
    doPlot(My,xl,ud.v,1)
end
if ~isempty(ud.callback)
    try eval(ud.callback);end
end
    

function doPatch(src,evt)
hf = gcf;
set(hf,'WindowButtonDownFcn',@doPatch,...
    'WindowButtonUpFcn',@UndoPatch,...
    'WindowButtonMotionFcn',{@movePatch},...
    'Pointer','fleur')

function doLb(src,evt)
hf = gcf;
set(hf,'WindowButtonDownFcn',@doLb,...
    'WindowButtonUpFcn',@UndoPatch,...
    'WindowButtonMotionFcn',{@moveLb},...
    'Pointer','left')

function doRb(src,evt)
hf = gcf;
set(hf,'WindowButtonDownFcn',@doRb,...
    'WindowButtonUpFcn',@UndoPatch,...
    'WindowButtonMotionFcn',{@moveRb},...
    'Pointer','right')

function UndoPatch(src,evt)
ha = gca;
hf = gcf;
set(hf,'WindowButtonMotionFcn',[],...
    'WindowButtonDownFcn',[],...
    'WindowButtonUpFcn',[],...
    'Pointer','arrow')
ud = get(ha,'userdata');
xw = get(ud.v.handles.axes,'xlim');
if ~ud.v.transpose
    My = ud.v.M*ud.y(xw(1):xw(2),:)';
else
    My = ud.v.M*ud.y(:,xw(1):xw(2));
end
doPlot(My,xw,ud.v,1);
if ~isempty(ud.callback)
    try eval(ud.callback);end
end


function movePatch(src,evt)
ha = gca;
cp = get(ha,'CurrentPoint');
ud = get(ha,'userdata');
xm = cp(1);
sw = diff(get(ud.v.handles.axes,'xlim'));
xl = xm + [-sw./2,+sw./2];
if xl(1) >= 1 && xl(2) <= ud.v.nt
    xl = round(xl);
    decim = max([1,round(sw./ud.v.ds)]);
    if ~ud.v.transpose
        My = ud.v.M*ud.y(xl(1):decim:xl(2),:)';
    else
        My = ud.v.M*ud.y(:,xl(1):decim:xl(2));
    end
    doPlot(My,xl,ud.v,decim)
elseif xl(2) > ud.v.nt
    xl = [ud.v.nt-sw,ud.v.nt];
    decim = max([1,round(sw./ud.v.ds)]);
    if ~ud.v.transpose
        My = ud.v.M*ud.y(xl(1):decim:xl(2),:)';
    else
        My = ud.v.M*ud.y(:,xl(1):decim:xl(2));
    end
    doPlot(My,xl,ud.v,decim)
elseif xl(1) < 1
    xl = [1,sw+1];
    decim = max([1,round(sw./ud.v.ds)]);
    if ~ud.v.transpose
        My = ud.v.M*ud.y(xl(1):decim:xl(2),:)';
    else
        My = ud.v.M*ud.y(:,xl(1):decim:xl(2));
    end
    doPlot(My,xl,ud.v,decim)
end

function moveLb(src,evt)
ha = gca;
cp = get(ha,'CurrentPoint');
cp = max([cp(1),1]);
ud = get(ha,'userdata');
xw = get(ud.v.handles.pa,'xdata');
if cp(1) <= xw(2) - ud.v.minSizeWindow ...
        && diff([cp(1) xw(2)]) <= ud.v.maxSizeWindow
    cp = [round(cp(1)) xw(2)];
    sw = diff(cp);
    decim = max([1,round(sw./ud.v.ds)]);
    if ~ud.v.transpose
        My = ud.v.M*ud.y(cp(1):decim:cp(2),:)';
    else
        My = ud.v.M*ud.y(:,cp(1):decim:cp(2));
    end
    doPlot(My,cp,ud.v,decim)
end

function moveRb(src,evt)
ha = gca;
cp = get(ha,'CurrentPoint');
ud = get(ha,'userdata');
cp = min([cp(1),ud.v.ind(end)]);
xw = get(ud.v.handles.pa,'xdata');
if xw(1) <= cp(1) - ud.v.minSizeWindow ...
        && diff([xw(1) cp(1)]) <= ud.v.maxSizeWindow
    cp = [xw(1) round(cp(1))];
    sw = diff(cp);
    decim = max([1,round(sw./ud.v.ds)]);
    if ~ud.v.transpose
        My = ud.v.M*ud.y(cp(1):decim:cp(2),:)';
    else
        My = ud.v.M*ud.y(:,cp(1):decim:cp(2));
    end
    doPlot(My,cp,ud.v,decim)
end

function doPlot(y,xw,v,decim)
for i=1:v.nc
    set(v.handles.hp(i),...
        'xdata',xw(1):decim:xw(2),...
        'ydata',y(i,:)+v.offset(i))
end
set(v.handles.axes,...
    'ylim',[v.mi v.ma],'xlim',xw);
set(v.handles.pa,...
    'xdata',[xw,fliplr(xw)]);
set(v.handles.lb,...
    'xdata',[xw(1) xw(1)]);
set(v.handles.rb,...
    'xdata',[xw(2) xw(2)]);
sw = diff(xw);
set(v.handles.hslider,...
    'value',mean(xw),...
    'min',max([1,sw/2-1]),...
    'max',max([v.nt,v.nt-sw/2+1]),...
    'sliderstep',.1*[sw/(v.nt-1) 4*sw/(v.nt-1)]);
