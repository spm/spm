function varargout = spm_mesh_render(action,varargin)
% Display a surface mesh & various utilities
% FORMAT H = spm_mesh_render('Disp',M,'PropertyName',propertyvalue)
% M        - a GIfTI filename/object or patch structure
% H        - structure containing handles of various objects
% Opens a new figure unless a 'parent' Property is provided with an axis
% handle.
%
% FORMAT H = spm_mesh_render(M)
% Shortcut to previous call format.
%
% FORMAT H = spm_mesh_render('ContextMenu',AX)
% AX       - axis handle or structure returned by spm_mesh_render('Disp',...)
%
% FORMAT H = spm_mesh_render('Overlay',AX,P)
% AX       - axis handle or structure given by spm_mesh_render('Disp',...)
% P        - data to be overlayed on mesh (see spm_mesh_project)
%
% FORMAT H = spm_mesh_render('ColourBar',AX,MODE)
% AX       - axis handle or structure returned by spm_mesh_render('Disp',...)
% MODE     - {['on'],'off'}
%
% FORMAT H = spm_mesh_render('ColourMap',AX,MAP)
% AX       - axis handle or structure returned by spm_mesh_render('Disp',...)
% MAP      - a colour map matrix
%
% FORMAT MAP = spm_mesh_render('ColourMap',AX)
% Retrieves the current colourmap.
%
% FORMAT spm_mesh_render('Register',AX,hReg)
% AX       - axis handle or structure returned by spm_mesh_render('Disp',...)
% hReg     - Handle of HandleGraphics object to build registry in.
% See spm_XYZreg for more information.
%__________________________________________________________________________
% Copyright (C) 2010-2011 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_render.m 6867 2016-09-12 15:04:44Z guillaume $


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
        if ischar(M) || isstruct(M), M = gifti(M); end
        if ~isfield(M,'vertices')
            try
                MM = M;
                M  = gifti(MM.private.metadata(1).value);
                try, M.cdata = MM.cdata(); end
            catch
                error('Cannot find a surface mesh to be displayed.');
            end
        end
        M = export(M,'patch');
        O = getOptions(varargin{2:end});
        
        %-Figure & Axis
        %------------------------------------------------------------------
        if isfield(O,'parent')
            H.axis   = O.parent;
            H.figure = ancestor(H.axis,'figure');
            figure(H.figure); axes(H.axis);
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
            'Clipping',         'off',...
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
        if isfield(M,'facevertexcdata')
            T = M.facevertexcdata;
        else
            T = [];
        end
        updateTexture(H,T);
        
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
        
        cmenu = uicontextmenu('Callback',{@myMenuCallback, H});
        
        uimenu(cmenu, 'Label','Inflate', 'Interruptible','off', ...
            'Callback',{@myInflate, H});
        
        uimenu(cmenu, 'Label','Overlay...', 'Interruptible','off', ...
            'Callback',{@myOverlay, H});
        
        uimenu(cmenu, 'Label','Image Sections...', 'Interruptible','off', ...
            'Callback',{@myImageSections, H});
        
        uimenu(cmenu, 'Label','Change geometry...', 'Interruptible','off', ...
            'Callback',{@myChangeGeometry, H});
        
        c = uimenu(cmenu, 'Label', 'Connected Components', 'Interruptible','off');
        C = getappdata(H.patch,'cclabel');
        for i=1:length(unique(C))
            uimenu(c, 'Label',sprintf('Component %d',i), 'Checked','on', ...
                'Callback',{@myCCLabel, H});
        end
        
        uimenu(cmenu, 'Label','Rotate', 'Checked','on', 'Separator','on', ...
            'Callback',{@mySwitchRotate, H});
        
        uimenu(cmenu, 'Label','Synchronise Views', 'Visible','off', ...
            'Checked','off', 'Tag','SynchroMenu', 'Callback',{@mySynchroniseViews, H});
        
        c = uimenu(cmenu, 'Label','View');
        uimenu(c, 'Label','Go to Y-Z view (right)',  'Callback', {@myView, H, [90 0]});
        uimenu(c, 'Label','Go to Y-Z view (left)',   'Callback', {@myView, H, [-90 0]});
        uimenu(c, 'Label','Go to X-Y view (top)',    'Callback', {@myView, H, [0 90]});
        uimenu(c, 'Label','Go to X-Y view (bottom)', 'Callback', {@myView, H, [-180 -90]});
        uimenu(c, 'Label','Go to X-Z view (front)',  'Callback', {@myView, H, [-180 0]});
        uimenu(c, 'Label','Go to X-Z view (back)',   'Callback', {@myView, H, [0 0]});
        
        uimenu(cmenu, 'Label','Colorbar', 'Callback', {@myColourbar, H});
        
        c = uimenu(cmenu, 'Label','Colormap');
        clrmp = {'hot' 'jet' 'gray' 'hsv' 'bone' 'copper' 'pink' 'white' ...
            'flag' 'lines' 'colorcube' 'prism' 'cool' 'autumn' ...
             'spring' 'winter' 'summer'};
        for i=1:numel(clrmp)
            uimenu(c, 'Label', clrmp{i}, 'Callback', {@myColourmap, H});
        end
        
        c = uimenu(cmenu, 'Label','Transparency');
        uimenu(c, 'Label','0%',  'Checked','on',  'Callback', {@myTransparency, H});
        uimenu(c, 'Label','20%', 'Checked','off', 'Callback', {@myTransparency, H});
        uimenu(c, 'Label','40%', 'Checked','off', 'Callback', {@myTransparency, H});
        uimenu(c, 'Label','60%', 'Checked','off', 'Callback', {@myTransparency, H});
        uimenu(c, 'Label','80%', 'Checked','off', 'Callback', {@myTransparency, H});
        
        uimenu(cmenu, 'Label','Data Cursor', 'Callback', {@myDataCursor, H});
        
        c = uimenu(cmenu, 'Label','Background Color');
        uimenu(c, 'Label','White',     'Callback', {@myBackgroundColor, H, [1 1 1]});
        uimenu(c, 'Label','Black',     'Callback', {@myBackgroundColor, H, [0 0 0]});
        uimenu(c, 'Label','Custom...', 'Callback', {@myBackgroundColor, H, []});
        
        uimenu(cmenu, 'Label','Save As...', 'Separator', 'on', ...
            'Callback', {@mySave, H});
        
        set(H.rotate3d,'enable','off');
        try, set(H.rotate3d,'uicontextmenu',cmenu); end
        try, set(H.patch,   'uicontextmenu',cmenu); end
        set(H.rotate3d,'enable','on');
        
        dcm_obj = datacursormode(H.figure);
        set(dcm_obj, 'Enable','off', 'SnapToDataVertex','on', ...
            'DisplayStyle','Window', 'Updatefcn',{@myDataCursorUpdate, H});
        
    %-Overlay
    %======================================================================
    case 'overlay'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if nargin < 3, varargin{2} = []; end
        updateTexture(H,varargin{2:end});
        
    %-Slices
    %======================================================================
    case 'slices'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if nargin < 3, varargin{2} = []; end
        renderSlices(H,varargin{2:end});
    
    %-ColourBar
    %======================================================================
    case {'colourbar', 'colorbar'}
        if isempty(varargin), varargin{1} = gca; end
        if length(varargin) == 1, varargin{2} = 'on'; end
        H   = getHandles(varargin{1});
        d   = getappdata(H.patch,'data');
        col = getappdata(H.patch,'colourmap');
        if strcmpi(varargin{2},'off')
            if isfield(H,'colourbar') && ishandle(H.colourbar)
                delete(H.colourbar);
                H = rmfield(H,'colourbar');
                setappdata(H.axis,'handles',H);
            end
            return;
        end
        if isempty(d) || ~any(d(:)), varargout = {H}; return; end
        if isempty(col), col = hot(256); end
        if ~isfield(H,'colourbar') || ~ishandle(H.colourbar)
            H.colourbar = colorbar('peer',H.axis);
            set(H.colourbar,'Tag','');
            set(get(H.colourbar,'Children'),'Tag','');
        end
        c(1:size(col,1),1,1:size(col,2)) = col;
        ic = findobj(H.colourbar,'Type','image');
        if size(d,1) > 1
            set(ic,'CData',c(1:size(d,1),:,:));
            set(ic,'YData',[1 size(d,1)]);
            set(H.colourbar,'YLim',[1 size(d,1)]);
            set(H.colourbar,'YTickLabel',[]);
        else
            set(ic,'CData',c);
            clim = getappdata(H.patch,'clim');
            if isempty(clim), clim = [false min(d) max(d)]; end
            set(ic,'YData',clim(2:3));
            set(H.colourbar,'YLim',clim(2:3));
        end
        setappdata(H.axis,'handles',H);
        
    %-ColourMap
    %======================================================================
    case {'colourmap', 'colormap'}
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1
            varargout = { getappdata(H.patch,'colourmap') };
            return;
        else
            setappdata(H.patch,'colourmap',varargin{2});
            d = getappdata(H.patch,'data');
            updateTexture(H,d);
        end
    
    %-CLim
    %======================================================================
    case 'clim'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1
            c = getappdata(H.patch,'clim');
            if ~isempty(c), c = c(2:3); end
            varargout = { c };
            return;
        else
            if isempty(varargin{2}) || any(~isfinite(varargin{2}))
                setappdata(H.patch,'clim',[false NaN NaN]);
            else
                setappdata(H.patch,'clim',[true varargin{2}]);
            end
            d = getappdata(H.patch,'data');
            updateTexture(H,d);
        end
        
    %-Register
    %======================================================================
    case 'register'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        hReg = varargin{2};
        xyz  = spm_XYZreg('GetCoords',hReg);
        hs   = myCrossBar('Create',H,xyz);
        set(hs,'UserData',hReg);
        spm_XYZreg('Add2Reg',hReg,hs,@myCrossBar);
        
    %-Debug
    %======================================================================
    case 'debug'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        set(H.patch,'EdgeColor','r');
        M = getappdata(H.patch,'patch');
        for i=1:size(M.vertices,1)
            text(M.vertices(i,1),M.vertices(i,2),M.vertices(i,3),...
                sprintf('%d',i));
        end
        
    %-Otherwise...
    %======================================================================
    otherwise
        try
            H = spm_mesh_render('Disp',action,varargin{:});
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
if ~nargin || isempty(H), H = gca; end
if ishandle(H) && ~isappdata(H,'handles')
    a = H; clear H;
    H.axis     = a;
    H.figure   = ancestor(H.axis,'figure');
    H.patch    = findobj(H.axis,'type','patch');
    H.light    = findobj(H.axis,'type','light');
    H.rotate3d = rotate3d(H.figure);
    setappdata(H.axis,'handles',H);
elseif ishandle(H)
    H = getappdata(H,'handles');
else
    H = getappdata(H.axis,'handles');
end

%==========================================================================
function myMenuCallback(obj,evt,H)
H = getHandles(H);

h = findobj(obj,'Label','Rotate');
if strcmpi(get(H.rotate3d,'Enable'),'on')
    set(h,'Checked','on');
else
    set(h,'Checked','off');
end

if numel(findobj('Tag','SPMMeshRender','Type','Patch')) > 1
    h = findobj(obj,'Tag','SynchroMenu');
    set(h,'Visible','on');
end

h = findobj(obj,'Label','Colorbar');
d = getappdata(H.patch,'data');
if isempty(d) || ~any(d(:)), set(h,'Enable','off'); else set(h,'Enable','on'); end
if isfield(H,'colourbar')
    if ishandle(H.colourbar)
        set(h,'Checked','on');
    else
        H = rmfield(H,'colourbar');
        set(h,'Checked','off');
    end
else
    set(h,'Checked','off');
end
setappdata(H.axis,'handles',H);

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
function varargout = myCrossBar(varargin)

switch lower(varargin{1})

    case 'create'
    %----------------------------------------------------------------------
    % hMe = myCrossBar('Create',H,xyz)
    H  = varargin{2};
    xyz = varargin{3};
    hold(H.axis,'on');
    hs = plot3(xyz(1),xyz(2),xyz(3),'Marker','+','MarkerSize',60,...
        'parent',H.axis,'Color',[1 1 1],'Tag','CrossBar','ButtonDownFcn',{});
    varargout = {hs};
    
    case 'setcoords'
    %----------------------------------------------------------------------
    % [xyz,d] = myCrossBar('SetCoords',xyz,hMe)
    hMe  = varargin{3};
    xyz  = varargin{2};
    set(hMe,'XData',xyz(1));
    set(hMe,'YData',xyz(2));
    set(hMe,'ZData',xyz(3));
    varargout = {xyz,[]};
    
    otherwise
    %----------------------------------------------------------------------
    error('Unknown action string')

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
function myTransparency(obj,evt,H)
t = 1 - sscanf(get(obj,'Label'),'%d%%') / 100;
set(H.patch,'FaceAlpha',t);
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');

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
function myColourbar(obj,evt,H)
y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')};
spm_mesh_render('Colourbar',H,toggle(get(obj,'Checked')));

%==========================================================================
function myColourmap(obj,evt,H)
spm_mesh_render('Colourmap',H,feval(get(obj,'Label'),256));

%==========================================================================
function mySynchroniseViews(obj,evt,H)
P = findobj('Tag','SPMMeshRender','Type','Patch');
v = get(H.axis,'cameraposition');
for i=1:numel(P)
    H = getappdata(ancestor(P(i),'axes'),'handles');
    set(H.axis,'cameraposition',v);
    axis(H.axis,'image');
    camlight(H.light);
end

%==========================================================================
function myDataCursor(obj,evt,H)
dcm_obj = datacursormode(H.figure);
set(dcm_obj, 'Enable','on', 'SnapToDataVertex','on', ...
    'DisplayStyle','Window', 'Updatefcn',{@myDataCursorUpdate, H});

%==========================================================================
function txt = myDataCursorUpdate(obj,evt,H)
pos = get(evt,'Position');
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Z: ',num2str(pos(3))]};
i = ismember(get(H.patch,'vertices'),pos,'rows');
txt = {['Node: ' num2str(find(i))] txt{:}};
d = getappdata(H.patch,'data');
if ~isempty(d) && any(d(:))
    if any(i), txt = {txt{:} ['T: ',num2str(d(i))]}; end
end
hMe = findobj(H.axis,'Tag','CrossBar');
if ~isempty(hMe)
    ws = warning('off');
    spm_XYZreg('SetCoords',pos,get(hMe,'UserData'));
    warning(ws);
end

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
    '*.idtf' 'IDTF files (*.idtf)';...
    '*.vtk' 'VTK files (*.vtk)'}, 'Save as');
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
        case {'.vtk','.vtp'}
            filterindex = 5;
        otherwise
            switch filterindex
                case 1
                    filename = [filename '.gii'];
                case 2
                    filename = [filename '.png'];
                case 3
                    filename = [filename '.dae'];
                case 4
                    filename = [filename '.idtf'];
                case 5
                    filename = [filename '.vtk'];
            end
    end
    switch filterindex
        case 1
            G = gifti(H.patch);
            [p,n,e] = fileparts(filename);
            [p,n,e] = fileparts(n);
            switch lower(e)
                case '.func'
                    save(gifti(getappdata(H.patch,'data')),...
                        fullfile(pathname, filename));
                case '.surf'
                    save(gifti(struct('vertices',G.vertices,'faces',G.faces)),...
                        fullfile(pathname, filename));
                case '.rgba'
                    save(gifti(G.cdata),fullfile(pathname, filename));
                otherwise
                    save(G,fullfile(pathname, filename));
            end
        case 2
            u  = get(H.axis,'units');
            set(H.axis,'units','pixels');
            p  = get(H.axis,'Position');
            r  = get(H.figure,'Renderer');
            hc = findobj(H.figure,'Tag','SPMMeshRenderBackground');
            if isempty(hc)
                c = get(H.figure,'Color');
            else
                c = get(hc,'Color');
            end
            h = figure('Position',p+[0 0 10 10], ...
                'InvertHardcopy','off', ...
                'Color',c, ...
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
            saveas(gifti(H.patch),fullfile(pathname, filename),'collada');
        case 4
            saveas(gifti(H.patch),fullfile(pathname, filename),'idtf');
        case 5
            saveas(gifti(H.patch),fullfile(pathname, filename),'vtk');
    end
end

%==========================================================================
function myDeleteFcn(obj,evt,renderer)
try, rotate3d(get(obj,'parent'),'off'); end
set(ancestor(obj,'figure'),'Renderer',renderer);

%==========================================================================
function myOverlay(obj,evt,H)
[P, sts] = spm_select(1,'\.img|\.nii|\.gii|\.mat','Select file to overlay');
if ~sts, return; end
spm_mesh_render('Overlay',H,P);

%==========================================================================
function myImageSections(obj,evt,H)
[P, sts] = spm_select(1,'image','Select image to render');
if ~sts, return; end
renderSlices(H,P);

%==========================================================================
function myChangeGeometry(obj,evt,H)
[P, sts] = spm_select(1,'mesh','Select new geometry mesh');
if ~sts, return; end
G = gifti(P);
if size(H.patch.Vertices,1) ~= size(G.vertices,1)
    error('Number of vertices must match.');
end
H.patch.Vertices = G.vertices;
H.patch.Faces = G.faces;

%==========================================================================
function renderSlices(H,P,pls)
if nargin <3
    pls = 0.05:0.2:0.9;
end
N   = nifti(P);
d   = size(N.dat);
pls = round(pls.*d(3));
hold(H.axis,'on');
for i=1:numel(pls)
    [x,y,z] = ndgrid(1:d(1),1:d(2),pls(i));
    f  = N.dat(:,:,pls(i));
    x1 = N.mat(1,1)*x + N.mat(1,2)*y + N.mat(1,3)*z + N.mat(1,4);
    y1 = N.mat(2,1)*x + N.mat(2,2)*y + N.mat(2,3)*z + N.mat(2,4);
    z1 = N.mat(3,1)*x + N.mat(3,2)*y + N.mat(3,3)*z + N.mat(3,4);
    surf(x1,y1,z1, repmat(f,[1 1 3]), 'EdgeColor','none', ...
        'Clipping','off', 'Parent',H.axis);
end
hold(H.axis,'off');
axis(H.axis,'image');

%==========================================================================
function C = updateTexture(H,v,col)

%-Get colourmap
%--------------------------------------------------------------------------
if nargin<3, col = getappdata(H.patch,'colourmap'); end
if isempty(col), col = hot(256); end
setappdata(H.patch,'colourmap',col);

%-Get curvature
%--------------------------------------------------------------------------
curv = getappdata(H.patch,'curvature');
if size(curv,2) == 1
    curv = 0.5 * repmat(curv,1,3) + 0.3 * repmat(~curv,1,3);
end
    
%-Project data onto surface mesh
%--------------------------------------------------------------------------
if nargin < 2, v = []; end
if ischar(v)
    [p,n,e] = fileparts(v);
    if strcmp([n e],'SPM.mat')
        swd = pwd;
        spm_figure('GetWin','Interactive');
        [SPM,v] = spm_getSPM(struct('swd',p));
        cd(swd);
    else
        try, spm_vol(v); catch, v = gifti(v); end;
    end
end
if isa(v,'gifti'), v = v.cdata; end
if isa(v,'file_array'), v = v(); end
if isempty(v)
    v = zeros(size(curv))';
elseif ischar(v) || iscellstr(v) || isstruct(v)
    v = spm_mesh_project(H.patch,v);
elseif isnumeric(v) || islogical(v)
    if size(v,2) == 1
        v = v';
    end
else
    error('Unknown data type.');
end
v(isinf(v)) = NaN;

setappdata(H.patch,'data',v);

%-Create RGB representation of data according to colourmap
%--------------------------------------------------------------------------
C = zeros(size(v,2),3);
clim = getappdata(H.patch, 'clim');
if isempty(clim), clim = [false NaN NaN]; end
mi = clim(2); ma = clim(3);
if any(v(:))
    if size(col,1)>3 && size(col,1) ~= size(v,1)
        if size(v,1) == 1
            if ~clim(1), mi = min(v(:)); ma = max(v(:)); end
            C = squeeze(ind2rgb(floor(((v(:)-mi)/(ma-mi))*size(col,1)),col));
        elseif isequal(size(v),[size(curv,1) 3])
            C = v; v = v';
        else
            if ~clim(1), mi = min(v(:)); ma = max(v(:)); end
            for i=1:size(v,1)
                C = C + squeeze(ind2rgb(floor(((v(i,:)-mi)/(ma-mi))*size(col,1)),col));
            end
        end
    else
        if ~clim(1), ma = max(v(:)); end
        for i=1:size(v,1)
            C = C + v(i,:)'/ma * col(i,:);
        end
    end
end
setappdata(H.patch, 'clim', [false mi ma]);

%-Build texture by merging curvature and data
%--------------------------------------------------------------------------
C = repmat(~any(v,1),3,1)' .* curv + repmat(any(v,1),3,1)' .* C;

set(H.patch, 'FaceVertexCData',C, 'FaceColor','interp');

%-Update the colourbar
%--------------------------------------------------------------------------
if isfield(H,'colourbar')
    spm_mesh_render('Colourbar',H);
end
