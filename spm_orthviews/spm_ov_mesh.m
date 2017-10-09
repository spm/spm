function ret = spm_ov_mesh(varargin)
% Mesh tool - plugin for spm_orthviews
%
% This routine is a plugin to spm_orthviews. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the MATLAB prompt.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% Torben Lund, Guillaume Flandin & Christian Gaser
% $Id: spm_ov_mesh.m 7183 2017-10-09 15:26:47Z guillaume $


switch lower(varargin{1})
    % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, ...
            'Label', 'Mesh',...
            'UserData',...
                struct('style','','width',1,'id',varargin{2}));
        item1 = uimenu(item0,...
            'Label','Select mesh(es)...',...
            'Tag','menu',...
            'Callback',@(hObj,event) mesh_menu(hObj,event,varargin{2}));
        item2 = uimenu(item0,...
            'Label','Options');
        item2_1 = uimenu(item2,...
            'Label','Line style...',...
            'UserData',struct('str','Line style', 'field','style'),...
            'Callback',@(hObj,event) mesh_options(hObj,event,item0));
        item2_2 = uimenu(item2,...
            'Label','Line width...',...
            'UserData',struct('str','Line width', 'field','width'),...
            'Callback',@(hObj,event) mesh_options(hObj,event,item0));
        ret = item0;
    case 'display'
        if nargin < 2
            varargin{2} = spm_input('Select image', '!+1', 'e');
        end
        mesh_add(varargin{2:end});
        mesh_delete(varargin{2});
        mesh_display(varargin{2});
    case 'redraw'
        mesh_redraw(varargin{2:end});
    otherwise
end


%==========================================================================
function mesh_options(hObj,event,hM)
Finter = spm_figure('FindWin', 'Interactive');
spm_input('!DeleteInputObj',Finter);
UDm = get(hObj,'UserData');
UDc = get(hM,'UserData');
if isnumeric(UDc.(UDm.field))
    N = Inf;
    T = 'e';
else
    N = Inf;
    T = 's';
end
if T == 's' && iscell(UDc.(UDm.field))
    UDc.(UDm.field) = '';
end
UDc.(UDm.field) = spm_input(UDm.str, '!+1', T, UDc.(UDm.field), N);
if T == 's'
    if ~isempty(strfind(UDc.(UDm.field),'{'))
        UDc.(UDm.field) = eval(UDc.(UDm.field));
    end
end
set(hM,'UserData',UDc);
try, mesh_redraw(UDc.id); end


%==========================================================================
function mesh_menu(hObj,event,i,varargin)

global st

mesh_add(i);
mesh_delete(i);
mesh_display(i);


%==========================================================================
function mesh_add(i,g)

global st

try
    m = st.vols{i}.mesh.meshes;
catch
    m = [];
end

if nargin < 2
    [g,sts] = spm_select([1 Inf],'mesh','Select mesh(es)...');
    if ~sts, st.vols{i}.mesh.meshes = m; return; end
end
g = gifti(g);
for j=1:numel(g)
    if ~isfield(g(j),'vertices') || ~isfield(g(j),'faces')
        error('Selected file must contain triangular meshes.');
    end
end

st.vols{i}.mesh.meshes = [m, g];

hM = findobj(st.vols{i}.ax{1}.cm,'Label','Mesh');
if numel(get(hM,'children')) < 3
    uimenu(hM,...
        'Label','Display 3D mesh(es)',...
        'Tag','tmp',...
        'Callback',@(hObj,event) mesh_render(i));
    uimenu(hM,...
        'Label','Remove mesh(es)',...
        'Tag','tmp',...
        'Callback',@(hObj,event) mesh_clear(i));
end


%==========================================================================
function mesh_clear(i)

global st

mesh_delete(i);
st.vols{i} = rmfield(st.vols{i},'mesh');
hM = findobj(st.vols{i}.ax{1}.cm,'Label','Mesh');
set(findobj(hM,'Tag','menu'),'Label','Select mesh(es)');
delete(findobj(get(hM,'Children'),'Tag','tmp'));


%==========================================================================
function mesh_display(i)

global st

if nargin < 1
    i = 1;
end
o = st.vols{i}.mesh.meshes;

try
    hM = findobj(st.vols{i}.ax{1}.cm,'Label','Mesh');
    set(findobj(hM,'Tag','menu'),'Label','Add mesh(es)');
    UD = get(hM,'UserData');
    linespec  = UD.style;
    linewidth = UD.width;
    set(hM,'UserData',UD);
end
if ~iscell(linespec)
    if ~isempty(linespec)
        disp('Wrong Linestyle format please use {''linestyle1'',''linestyle2'',...}.');
        UD.style = '';
        set(hM,'UserData',UD);
    end
    
    linespec = {'b-' 'g-' 'r-' 'c-' 'm-' 'y-' 'k-' 'w-' ...
        'b-.' 'g-.' 'r-.' 'c-.' 'm-.' 'y-.' 'k-.' 'w-.' ...
        'b--' 'g--' 'r--' 'c--' 'm--' 'y--' 'k--' 'w--' };
    if numel(o) <= 24
        linespec = linespec(1:numel(o));
    else
        % If there are more than 24 surfaces plot them all in red
        linespec = repmat({'r-'},1,numel(o));
    end
end
if numel(linespec) < numel(o)
    disp('Meshes with unspecified line style are displayed in black.')
    linespec(end+1:numel(o)) = {'k-'};
    UD.style = linespec;
    set(hM,'UserData',UD);
end
if numel(linewidth) < numel(o)
    linewidth(end+1:numel(o)) = linewidth(end);
    UD.width  = linewidth;
    set(hM,'UserData',UD);
end

% This requires more work: does not work in voxel space, ie handling of
% st.Space is incomplete
bb = st.bb;
is = inv(st.Space);
lh = {};
for n=1:numel(o)
    vert = (st.vols{i}.premul(1:3,:)*[o(n).vertices';ones(1,size(o(n).vertices,1))])';
    % 1
    set(st.vols{i}.ax{1}.ax,'NextPlot','add');
    s = spm_mesh_contour(...
        struct('faces',o(n).faces,'vertices',vert),...
        st.centre(3));
    for j=1:numel(s)
        c = is*[s(j).xdata; s(j).ydata;ones(2,size(s(j).xdata,2))];
        c(2,:) = c(2,:)-bb(1,2)+1;
        c(1,:) = c(1,:)-bb(1,1)+1;
        lh{end+1} = plot(st.vols{i}.ax{1}.ax,c(1,:),c(2,:),...
            linespec{n},'LineWidth',linewidth(n));
    end
    % 2
    set(st.vols{i}.ax{2}.ax,'NextPlot','add');
    s = spm_mesh_contour(...
        struct('faces',o(n).faces,'vertices',vert(:,[3 1 2])),...
        st.centre(2));
    for j=1:numel(s)
        c = is*[s(j).xdata;ones(1,size(s(j).xdata,2));s(j).ydata;ones(1,size(s(j).ydata,2))];
        c(3,:) = c(3,:)-bb(1,1)+1;
        c(1,:) = c(1,:)-bb(1,3)+1;
        lh{end+1} = plot(st.vols{i}.ax{2}.ax,c(3,:),c(1,:),...
            linespec{n},'LineWidth',linewidth(n));
    end
    % 3
    if st.mode == 0
        warning('Not handled.');
    else
        set(st.vols{i}.ax{3}.ax,'NextPlot','add');
        s = spm_mesh_contour(...
            struct('faces',o(n).faces,'vertices',vert(:,[2 3 1])),...
            st.centre(1));
        for j=1:numel(s)
            c = is*[ones(1,size(s(j).xdata,2));[s(j).xdata; s(j).ydata];ones(1,size(s(j).ydata,2))];
            c(2,:) = -c(2,:)+bb(2,2)+1;
            c(3,:) = c(3,:)-bb(1,3)+1;
            lh{end+1} = plot(st.vols{i}.ax{3}.ax,c(2,:),c(3,:),...
                linespec{n},'LineWidth',linewidth(n));
        end
    end
end

set(cat(1,lh{:}),'HitTest','off');

st.vols{i}.mesh.handles = lh;


%==========================================================================
function mesh_redraw(i,varargin) %i, TM0, TD, CM0, CD, SM0, SD

mesh_delete(i);
mesh_display(i);


%==========================================================================
function mesh_delete(i)

global st

try, delete(cat(1,st.vols{i}.mesh.handles{:})); end


%==========================================================================
function mesh_render(i)

global st

m = st.vols{i}.mesh.meshes;
for j=1:numel(m)
    spm_mesh_render(m(j));
end
% data cursor should be linked to orthviews through spm_XYZreg
% see spm_ovhelper_3Dreg.m
