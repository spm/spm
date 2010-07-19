function varargout = spm_mesh_render(action,varargin)
% Display a surface mesh & various utilities
% FORMAT H = spm_mesh_render('Disp',M)
% M        - a GIfTI filename/object or patch structure
% H        - structure containing handles of various objects
%
% FORMAT H = spm_mesh_render(M)
% Shortcut to previous call format
%
% FORMAT H = spm_mesh_render('ContextMenu',AX)
% AX       - axis handle or structure given by spm_mesh_render('Disp',...)
%
% FORMAT H = spm_mesh_render('AddOverlay',AX,P)
% AX       - axis handle or structure given by spm_mesh_render('Disp',...)
% P        - data to be overlayed on mesh (see spm_mesh_project)
%
% FORMAT H = spm_mesh_render('AddColourBar',AX)
% AX       - axis handle or structure given by spm_mesh_render('Disp',...)
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_render.m 4003 2010-07-19 18:22:38Z guillaume $


%-Input parameters
%--------------------------------------------------------------------------
if ~nargin, action = 'Disp'; end

if ~ischar(action)
    varargin = {action varargin{:}};
    action   = 'Disp';
end

varargout = {[]};

%-Action
%--------------------------------------------------------------------------
switch lower(action)
    
    %-Display
    %======================================================================
    case 'disp'
        if isempty(varargin)
            [M, sts] = spm_select(1,'mesh','Select surface mesh file');
            if ~sts, return; end
        else
            M = varargin{1};
        end
        if ischar(M), M = export(gifti(M),'patch'); end
        O = getOptions(varargin{2:end});
        
        %-Figure & Axis
        %------------------------------------------------------------------
        if isfield(O,'parent')
            H.axis   = O.parent;
            H.figure = ancestor(H.axis,'figure');
        else
            H.figure = figure('Color',[1 1 1]);
            H.axis   = axes('Parent',H.figure);
            set(H.axis,'Visible','off');
        end
        renderer = get(H.figure,'Renderer');
        set(H.figure,'Renderer','OpenGL');
        
        %-Patch
        %------------------------------------------------------------------
        P = struct('vertices',M.vertices, 'faces',M.faces);
        H.patch = patch(P,...
            'FaceColor',        [0.6 0.6 0.6],...
            'EdgeColor',        'none',...
            'FaceLighting',     'phong',...
            'SpecularStrength', 0.7,...
            'AmbientStrength',  0.1,...
            'DiffuseStrength',  0.7,...
            'SpecularExponent', 10,...
            'DeleteFcn',        {@myDeleteFcn, renderer},...
            'Visible',          'off',...
            'Tag',              'SPMMeshRender',...
            'Parent',           H.axis);
        setappdata(H.patch,'patch',P);
        
        %-Label connected components of the mesh
        %------------------------------------------------------------------
        C = spm_mesh_label(P);
        setappdata(H.patch,'cclabel',C);
        
        %-Compute mesh curvature
        %------------------------------------------------------------------
        curv = spm_mesh_curvature(P) > 0;
        setappdata(H.patch,'curvature',curv);
        
        %-Apply texture to mesh
        %------------------------------------------------------------------
        updateTexture(H,[]);
        
        %-Set viewpoint, light and manipulation options
        %------------------------------------------------------------------
        axis(H.axis,'image');
        axis(H.axis,'off');
        view(H.axis,[-90 0]);
        material(H.figure,'dull');
        H.light = camlight; set(H.light,'Parent',H.axis);
        
        H.rotate3d = rotate3d(H.axis);
        set(H.rotate3d,'Enable','on');
        set(H.rotate3d,'ActionPostCallback',{@myPostCallback, H});
        %try
        %    setAllowAxesRotate(H.rotate3d, ...
        %        setxor(findobj(H.figure,'Type','axes'),H.axis), false);
        %end
        
        %-Store handles
        %------------------------------------------------------------------
        setappdata(H.axis,'handles',H);
        set(H.patch,'Visible','on');
        
        %-Add context menu
        %------------------------------------------------------------------
        spm_mesh_render('ContextMenu',H);
        
    %-Context Menu
    %======================================================================
    case 'contextmenu'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if ~isempty(get(H.patch,'UIContextMenu')), return; end
        
        cmenu = uicontextmenu('Callback',{@myMenuCallback});
        
        c1 = uimenu(cmenu, 'Label','Inflate', 'Interruptible','off', ...
            'Callback',{@myInflate, H});
        
        c1 = uimenu(cmenu, 'Label','Overlay...', 'Interruptible','off', ...
            'Callback',{@myOverlay, H});
        
        c2 = uimenu(cmenu, 'Label', 'Connected Components', 'Interruptible','off');
        C=getappdata(H.patch,'cclabel');
        for i=1:length(unique(C))
            uimenu(c2, 'Label',sprintf('Component %d',i), 'Checked','on', ...
                'Callback',{@myCCLabel, H});
        end
        
        c3 = uimenu(cmenu, 'Label','Rotate', 'Checked','on', 'Separator','on', ...
            'Callback',{@mySwitchRotate, H});
        
        c3 = uimenu(cmenu, 'Label','Synchronise Views', 'Visible','off', ...
            'Checked','off', 'Tag','SynchroMenu', 'Callback',{@mySynchroniseViews, H});
        
        c4 = uimenu(cmenu, 'Label','View');
        uimenu(c4, 'Label','Go to Y-Z view (right)',  'Callback', {@myView, H, [90 0]});
        uimenu(c4, 'Label','Go to Y-Z view (left)',   'Callback', {@myView, H, [-90 0]});
        uimenu(c4, 'Label','Go to X-Y view (top)',    'Callback', {@myView, H, [0 90]});
        uimenu(c4, 'Label','Go to X-Y view (bottom)', 'Callback', {@myView, H, [-180 -90]});
        uimenu(c4, 'Label','Go to X-Z view (front)',  'Callback', {@myView, H, [-180 0]});
        uimenu(c4, 'Label','Go to X-Z view (back)',   'Callback', {@myView, H, [0 0]});
        
        c5 = uimenu(cmenu, 'Label','Transparency');
        uimenu(c5, 'Label','0%',  'Checked','on',  'Callback', {@myTransparency, H});
        uimenu(c5, 'Label','20%', 'Checked','off', 'Callback', {@myTransparency, H});
        uimenu(c5, 'Label','40%', 'Checked','off', 'Callback', {@myTransparency, H});
        uimenu(c5, 'Label','60%', 'Checked','off', 'Callback', {@myTransparency, H});
        uimenu(c5, 'Label','80%', 'Checked','off', 'Callback', {@myTransparency, H});
        
        c6 = uimenu(cmenu, 'Label','Background Color');
        uimenu(c6, 'Label','White',     'Callback', {@myBackgroundColor, H, [1 1 1]});
        uimenu(c6, 'Label','Black',     'Callback', {@myBackgroundColor, H, [0 0 0]});
        uimenu(c6, 'Label','Custom...', 'Callback', {@myBackgroundColor, H, []});
        
        c7 = uimenu(cmenu, 'Label','Save As...', 'Separator', 'on', ...
            'Callback', {@mySave, H});
        
        try, set(H.rotate3d,'uicontextmenu',cmenu); end
        try, set(H.patch,   'uicontextmenu',cmenu); end
        
    %-AddOverlay
    %======================================================================
    case 'addoverlay'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if nargin < 3, varargin{2} = []; end
        if nargin < 4, varargin{3} = []; end
        updateTexture(H,varargin{2},varargin{3});
        
    %-AddColourBar
    %======================================================================
    case 'addcolourbar'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        d = getappdata(H.patch,'data');
        col = getappdata(H.patch,'colourmap');
        if isempty(d) || ~any(d(:)), return; end
        if isempty(col), col = hot(256); end
        H.colourbar = colorbar('peer',H.axis);
        c(1:size(col,1),1,1:size(col,2)) = col;
        set(get(H.colourbar,'child'),'CData',c);
        set(get(H.colourbar,'child'),'YData',[min(d) max(d)]);
        set(H.colourbar,'YLim',[min(d) max(d)]);
        set(H.colourbar,'Tag','');
        set(get(H.colourbar,'child'),'Tag','');
        setappdata(H.axis,'handles',H);
        
    %-Otherwise...
    %======================================================================
    otherwise
        try
            spm_mesh_render('Disp',action,varargin{:});
        catch
            error('Unknown action.');
        end
end

varargout = {H};


%==========================================================================
function O = getOptions(varargin)
O = [];
if ~nargin
    return;
elseif nargin == 1 && isstruct(varargin{1})
    for i=fieldnames(varargin{1})
        O.(lower(i{1})) = varargin{1}.(i{1});
    end
elseif mod(nargin,2) == 0
    for i=1:2:numel(varargin)
        O.(lower(varargin{i})) = varargin{i+1};
    end
else
    error('Invalid list of property/value pairs.');
end

%==========================================================================
function H = getHandles(H)
if ishandle(H) && ~isappdata(H,'handles')
    H.axis     = H;
    H.figure   = ancestor(H.axis,'figure');
    H.patch    = findobj(H.axis,'type','patch');
    H.light    = findobj(H.axis,'type','light');
    H.rotate3d = rotate3d(H.figure);
    setappdata(H.axis,'handles',H);
elseif ishandle(H)
    H = getappdata(H,'handles');
else
    %H = H;
end

%==========================================================================
function myPostCallback(obj,evt,H)
P = findobj('Tag','SPMMeshRender','Type','Patch');
if numel(P) == 1
    camlight(H.light);
else
    for i=1:numel(P)
        H = getappdata(ancestor(P(i),'axes'),'handles');
        camlight(H.light);
    end
end
    
%==========================================================================
function myMenuCallback(obj,evt)
if numel(findobj('Tag','SPMMeshRender','Type','Patch')) > 1
    h = findobj(obj,'Tag','SynchroMenu');
    if ~isempty(h)
        set(h,'Visible','on');
    end
end

%==========================================================================
function myInflate(obj,evt,H)
spm_mesh_inflate(H.patch,Inf,1);
axis(H.axis,'image');

%==========================================================================
function myCCLabel(obj,evt,H)
C   = getappdata(H.patch,'cclabel');
F   = get(H.patch,'Faces');
ind = sscanf(get(obj,'Label'),'Component %d');
V   = get(H.patch,'FaceVertexAlphaData');
Fa  = get(H.patch,'FaceAlpha');
if ~isnumeric(Fa)
    if ~isempty(V), Fa = max(V); else Fa = 1; end
    if Fa == 0, Fa = 1; end
end
if isempty(V) || numel(V) == 1
    Ve = get(H.patch,'Vertices');
    if isempty(V) || V == 1
        V = Fa * ones(size(Ve,1),1);
    else
        V = zeros(size(Ve,1),1);
    end
end
if strcmpi(get(obj,'Checked'),'on')
    V(reshape(F(C==ind,:),[],1)) = 0;
    set(obj,'Checked','off');
else
    V(reshape(F(C==ind,:),[],1)) = Fa;
    set(obj,'Checked','on');
end
set(H.patch, 'FaceVertexAlphaData', V);
if all(V)
    set(H.patch, 'FaceAlpha', Fa);
else
    set(H.patch, 'FaceAlpha', 'interp');
end

%==========================================================================
function mySwitchRotate(obj,evt,H)
if strcmpi(get(H.rotate3d,'enable'),'on')
    set(H.rotate3d,'enable','off');
    set(obj,'Checked','off');
else
    set(H.rotate3d,'enable','on');
    set(obj,'Checked','on');
end

%==========================================================================
function myView(obj,evt,H,varargin)
view(H.axis,varargin{1});
axis(H.axis,'image');
camlight(H.light);

%==========================================================================
function mySynchroniseViews(obj,evt,H,varargin)
P = findobj('Tag','SPMMeshRender','Type','Patch');
v = get(H.axis,'cameraposition');
for i=1:numel(P)
    H = getappdata(ancestor(P(i),'axes'),'handles');
    set(H.axis,'cameraposition',v);
    axis(H.axis,'image');
    camlight(H.light);
end

%==========================================================================
function myTransparency(obj,evt,H)
t = 1 - sscanf(get(obj,'Label'),'%d%%') / 100;
set(H.patch,'FaceAlpha',t);
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');

%==========================================================================
function myBackgroundColor(obj,evt,H,varargin)
if isempty(varargin{1})
    c = uisetcolor(H.figure, ...
        'Pick a background color...');
    if numel(c) == 1, return; end
else
    c = varargin{1};
end
h = findobj(H.figure,'Tag','SPMMeshRenderBackground');
if isempty(h)
    set(H.figure,'Color',c);
else
    set(h,'Color',c);
end

%==========================================================================
function mySave(obj,evt,H)
[filename, pathname, filterindex] = uiputfile({...
    '*.gii' 'GIfTI files (*.gii)'; ...
    '*.png' 'PNG files (*.png)';...
    '*.dae' 'Collada files (*.dae)';...
    '*.idtf' 'IDTF files (*.idtf)'}, 'Save as');
if ~isequal(filename,0) && ~isequal(pathname,0)
    [pth,nam,ext] = fileparts(filename);
    switch ext
        case '.gii'
            filterindex = 1;
        case '.png'
            filterindex = 2;
        case '.dae'
            filterindex = 3;
        case '.idtf'
            filterindex = 4;
        otherwise
            switch filterindex
                case 1
                    filename = [filename '.gii'];
                case 2
                    filename = [filename '.png'];
                case 3
                    filename = [filename '.dae'];
            end
    end
    switch filterindex
        case 1
            save(gifti(H.patch),fullfile(pathname, filename));
        case 2
            u = get(H.axis,'units');
            set(H.axis,'units','pixels');
            p = get(H.axis,'Position');
            r = get(getappdata(obj,'fig'),'Renderer');
            h = figure('Position',p+[0 0 10 10], ...
                'InvertHardcopy','off', ...
                'Color',get(H.figure,'Color'), ...
                'Renderer',r);
            copyobj(H.axis,h);
            set(H.axis,'units',u);
            set(get(h,'children'),'visible','off');
            %a = get(h,'children');
            %set(a,'Position',get(a,'Position').*[0 0 1 1]+[10 10 0 0]);       
            if isdeployed
                deployprint(h, '-dpng', '-opengl', fullfile(pathname, filename));
            else
                print(h, '-dpng', '-opengl', fullfile(pathname, filename));
            end
            close(h);
            set(getappdata(obj,'fig'),'renderer',r);
        case 3
            save(gifti(H.patch),fullfile(pathname, filename),'collada');
        case 4
            save(gifti(H.patch),fullfile(pathname, filename),'idtf');
    end
end

%==========================================================================
function myDeleteFcn(obj,evt,renderer)
try, rotate3d(get(obj,'parent'),'off'); end
set(ancestor(obj,'figure'),'Renderer',renderer);

%==========================================================================
function myOverlay(obj,evt,H)
[P, sts] = spm_select(1,'image','Select image to overlay');
if ~sts, return; end
spm_mesh_render('AddOverlay',H,P);

%==========================================================================
function C = updateTexture(H,v,col)

if nargin <3 || isempty(col), col = hot(256); end
setappdata(H.patch,'colourmap',col);

%-Get curvature
%--------------------------------------------------------------------------
curv = getappdata(H.patch,'curvature');
if size(curv,2) == 1
    curv = 0.5 * repmat(curv,1,3) + 0.3 * repmat(~curv,1,3);
end
    
%-Project data onto surface mesh
%--------------------------------------------------------------------------
if isstruct(v) && isfield(v,'FWHM') %-xSPM structure
    v = struct('XYZ',v.XYZ, 't',v.Z, 'mat',v.M, 'dim',v.DIM);
end
if isempty(v)
    v = zeros(size(curv))';
elseif ischar(v) || iscellstr(v) || isstruct(v)
    v = spm_mesh_project(H.patch,v);
end

setappdata(H.patch,'data',v);

%-Build texture by merging curvature and data
%--------------------------------------------------------------------------
C = zeros(size(v,2),3);
if any(v(:))
    if size(col,1)>3
        if size(v,1) == 1
            C = squeeze(ind2rgb(floor(v(:)/max(v(:))*size(col,1)),col));
        else
            C = v; v = v';
        end
    else
        m = max(v(:));
        for i=1:size(v,1)
            C = C + v(i,:)'/m * col(i,:);
        end
    end
end
C = repmat(~any(v,1),3,1)' .* curv + repmat(any(v,1),3,1)' .* C;

set(H.patch, 'FaceVertexCData',C, 'FaceColor','interp');

% if isfield(M,'cdata') && ~isempty(M.cdata)
%     if any(size(M.cdata)==1), M.cdata = M.cdata(:); end
%     if size(M.cdata,1)==3,    M.cdata = M.cdata';   end
%     set(H.patch, 'FaceVertexCData',M.cdata, 'FaceColor','interp');
% else
%     set(H.patch, 'FaceVertexCData',curv, 'FaceColor','interp');
% end
