function ret = spm_ov_contour(varargin)
% Contour tool - plugin for spm_orthviews
%
% This routine is a plugin to spm_orthviews. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the MATLAB prompt.
%__________________________________________________________________________
% Copyright (C) 2014-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_ov_contour.m 6774 2016-04-20 12:25:08Z guillaume $


switch lower(varargin{1})
    % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, ...
            'Label', 'Contour',...
            'UserData',struct('nblines',3,'style','r-','id',varargin{2}));
        item1 = uimenu(item0,...
            'Label','Display onto');
        item1_1 = uimenu(item1,...
            'Label','all',...
            'Callback',@(hObj,event) contour_menu(hObj,event,varargin{2},Inf));
        item1_2 = uimenu(item1,...
            'Label','all but current',...
            'Callback',@(hObj,event) contour_menu(hObj,event,varargin{2},NaN));
        item1_3 = uimenu(item1,...
            'Label','select...',...
            'Callback',@(hObj,event) contour_menu(hObj,event,varargin{2}));
        item2 = uimenu(item0,...
            'Label','Options');
        item2_1 = uimenu(item2,...
            'Label', 'Number of lines...',...
            'UserData',struct('str','Number of lines', 'field','nblines'),...
            'Callback',@(hObj,event) contour_options(hObj,event,item0));
        item2_2 = uimenu(item2,...
            'Label','Line style...',...
            'UserData',struct('str','Line style', 'field','style'),...
            'Callback',@(hObj,event) contour_options(hObj,event,item0));
        %'Callback',@(hObj,event) contour_display(hObj,event,varargin{2}));
        ret = item0;
    case 'display'
        if nargin == 1
            varargin{2} = spm_input('Select image', '!+1', 'e');
        end
        contour_display(varargin{2:end});
    case 'redraw'
        contour_redraw(varargin{2:end});
    otherwise
end


%==========================================================================
function contour_options(hObj,event,hM)
Finter = spm_figure('FindWin', 'Interactive');
spm_input('!DeleteInputObj',Finter);
UDm = get(hObj,'UserData');
UDc = get(hM,'UserData');
if isnumeric(UDc.(UDm.field))
    N = 1;
    T = 'e';
else
    N = Inf;
    T = 's';
end
UDc.(UDm.field) = spm_input(UDm.str, '!+1', T, UDc.(UDm.field), N);
set(hM,'UserData',UDc);
try, contour_redraw(UDc.id); end


%==========================================================================
function contour_menu(hObj,event,i,varargin)

global st

allOpts = get(get(hObj,'parent'),'children');
current = strcmp(get(hObj,'Checked'),'on');
if any(strcmp(get(allOpts,'Checked'),'on'))
    set(allOpts, 'Checked','off');
    contour_delete(i);
    st.vols{i} = rmfield(st.vols{i},'contour');
end
if ~current
    set(hObj, 'Checked','on');
    contour_display(i,varargin{:}); 
end


%==========================================================================
function contour_display(i,o)

global st

if nargin < 2
    o = spm_input('Select image(s)', '!+1', 'e', ...
        num2str(spm_orthviews('valid_handles')));
    o = intersect(spm_orthviews('valid_handles'),o);
elseif isinf(o)
    o = spm_orthviews('valid_handles');
elseif isnan(o)
    o = setxor(spm_orthviews('valid_handles'),i);
end

try
    hM = findobj(st.vols{i}.ax{1}.cm,'Label','Contour');
    UD = get(hM,'UserData');
    nblines = UD.nblines;
    linestyle = UD.style;
catch
    nblines = 3;
    linestyle = 'r-';
end
linewidth = 1;

contour_delete(i);

lh = {};
sw = warning('off','MATLAB:contour:ConstantData');
for d = 1:3
    CData = sqrt(sum(get(st.vols{i}.ax{d}.d,'CData').^2, 3));
    CData(isinf(CData)) = NaN;
    CData(isnan(CData)) = 0;
    for h = o(:)'
        set(st.vols{h}.ax{d}.ax,'NextPlot','add');
        [C,lh{end+1}] = ...
            contour(st.vols{h}.ax{d}.ax,CData,...
            nblines,linestyle,'LineWidth',linewidth);
    end
end
warning(sw);
set(cat(1,lh{:}),'HitTest','off');

st.vols{i}.contour.images = o;
st.vols{i}.contour.handles = lh;


%==========================================================================
function contour_redraw(i,varargin) %i, TM0, TD, CM0, CD, SM0, SD

global st

contour_delete(i);
contour_display(i,st.vols{i}.contour.images);


%==========================================================================
function contour_delete(i)

global st

try, delete(cat(1,st.vols{i}.contour.handles{:})); end
