function [out] = spm_eeg_render(m,options)
% Visualisation routine for the cortical surface
% FORMAT [out] = spm_eeg_render(m,options)
%
% INPUT:
% m       - patch structure (with fields .faces et .vertices)
%           or GIFTI format filename
% options - structure with optional fields:
%           .texture      - texture to be projected onto the mesh
%           .clusters     - cortical parcellation (cell variable containing
%                           the vertex indices of each cluster) 
%           .clustersName - name of the clusters
%           .figname      - name to be given to the figure
%           .ParentAxes   - handle of the axes within which the mesh should
%                           be displayed
%           .hfig         - handle of existing figure. If this option is
%                           provided, then spm_eeg_render adds the (textured)
%                           mesh to the figure hfig, and a control for its
%                           transparency.
%
% OUTPUT:
% out     - structure with fields:
%           .hfra    - frame structure for movie building
%           .handles - structure containing the handles of the created
%                      uicontrols and mesh objects
%           .m       - the structure used to create the mesh.
%__________________________________________________________________________
%
% This function is a visualisation routine, mainly for texture and
% clustering on the cortical surface.
% NB: The texture and the clusters cannot be visualised at the same time.
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_render.m 6808 2016-06-13 16:48:30Z guillaume $


%----------------------------------------------------------------------%
%------------- Common features for any visualisation ------------------%
%----------------------------------------------------------------------%

% Check mesh format
try
    if ~isstruct(m)
        m = export(gifti(m),'patch');
    end
catch
    if ischar(m), fprintf('File: %s\n',m); end
    error('Unknown mesh format.')
end


% Default options
handles.fi = figure(...
    'visible','off',...
    'color',ones(1,3),...
    'NumberTitle','Off',...
    'Name','Mesh visualisation',...
    'tag','spm_eeg_render');
ns = 0;
texture = 'none';
clusters = 'none';
subplotBIN = 0;
addMesh = 0;
tag = '';
visible = 'on';
ParentAxes = axes('parent',handles.fi);
try, options; catch, options = [];end
% Now get options
if ~isempty(options)
    % get texture if provided
    try texture = options.texture;end
    % get ParentAxes
    try ParentAxes = options.ParentAxes;end
    % get tag
    try tag = options.tag;end
    % get flag for visibility: useful for displaying all objects at once
    try visible = options.visible;end
    % get custers if provided
    try
        clusters = options.clusters;
        IND = zeros(1,length(m.vertices));
        K = length(clusters);
        for k = 1:K
            IND(clusters{k}) = k+1./K;
        end
        texture = IND';
    end
    % get figname if provided
    try
        set(handles.fi,'NumberTitle','Off','Name',options.figname);
    end
    % get figure handle (should be parent of ParentAxes)
    try
        figure(options.hfig)
        if isempty(ParentAxes)
            ParentAxes = axes('parent',options.hfig,...
                'nextplot','add');
        end
        close(handles.fi);
        handles.fi = options.hfig;
        addMesh = 1;
        try % get number of transparency sliders in current figure...
            hh=get(handles.fi,'children');
            ns=length(findobj(hh,'userdata','tag_UIC_transparency'));
        catch
            ns=1;
        end
    end
    
end

handles.ParentAxes = ParentAxes;
oldRenderer = get(handles.fi,'renderer');
try
    if ismac
        set(handles.fi,'renderer','zbuffer');
    else
        set(handles.fi,'renderer','OpenGL');
    end
catch
    set(handles.fi,'renderer','OpenGL');
end

% Plot mesh and texture/clusters
if isequal(texture,'none')
    figure(handles.fi)
    handles.p = patch(m,...
        'facecolor', [.5 .5 .5], 'EdgeColor', 'none',...
        'FaceLighting','gouraud',...
        'parent',ParentAxes,...
        'userdata',oldRenderer,...
        'visible',visible,...
        'tag',tag);
else
    texture = texture(:);
    figure(handles.fi)
    if isequal(length(texture),length(m.vertices))
        handles.p = patch(m,...
            'facevertexcdata',texture,...
            'facecolor','interp',...
            'EdgeColor', 'none',...
            'FaceLighting','gouraud',...
            'parent',ParentAxes,...
            'userdata',oldRenderer,...
            'visible',visible,...
            'tag',tag,...
            'deleteFcn',@doDelMesh);
        col = colormap(ParentAxes,jet(256));
        udd.tex = texture;
        udd.cax = caxis(ParentAxes);
    else
        texture = 'none';
        disp('Warning: size of texture does not match number of vertices!')
        handles.p = patch(m,'facecolor', [.5 .5 .5], 'EdgeColor', 'none',...
            'parent',ParentAxes,...
            'userdata',oldRenderer,...
            'visible',visible,...
            'tag',tag,...
            'deleteFcn',@doDelMesh);
    end
end


daspect(ParentAxes,[1 1 1]);
axis(ParentAxes,'tight');
axis(ParentAxes,'off')
try,camva(ParentAxes,'auto');end
set(ParentAxes,'view',[25,45]);

% build internal userdata structure
udd.p = handles.p;


%----------------------------------------------------------------------%
%---------------------- GUI tools and buttons -------------------------%
%----------------------------------------------------------------------%


% Transparancy sliders
pos = [20 100 20 245];
pos(1) = pos(1) + ns.*25;
handles.transp = uicontrol(handles.fi,...
    'style','slider',...
    'position',pos,...
    'min',0,...
    'max',1,...
    'value',1,...
    'sliderstep',[0.01 0.05],...
    'userdata',handles.p,...
    'tooltipstring',['mesh #',num2str(ns+1),' transparency control'],...
    'callback',{@doTransp},...
    'BusyAction','cancel',...
    'Interruptible','off',...
    'visible',visible,...
    'tag',tag);
set(handles.transp,'units','normalized')
handles.tag = uicontrol(handles.fi,...
    'style','text',...
    'visible','off',...
    'tag',tag,...
    'userdata','tag_UIC_transparency');

udd.transp = handles.transp;


% Clustering buttons and popup menu
if ~isequal(clusters,'none')
    if subplotBIN
        subplot(2,1,1)
    end
%     set(p,'FaceColor','flat');
    col=lines;
    nc = floor(256./K);
    col = [repmat([0.8157 0.6666 0.5762],nc/2,1);...
        kron(col(1:K,:),ones(nc,1))];
    if K > 1
        col(end-nc/2:end,:) = [];
    end
    colormap(ParentAxes,col);
    tex = zeros(length(m.vertices),length(clusters)+1);
    tex(:,1) = texture;
    string = cell(length(clusters)+1,1);
    string{1} = 'all clusters';
    for i = 1:length(clusters)
        if ~isfield(options,'clustersName')
            string{i+1} = ['cluster ',num2str(i)];
        else
            string{i+1} = options.clustersName{i};
        end
        tex(clusters{i},i+1) = 1;
    end
    udd.tex = tex;
    udd.tex0 = tex;
    udd.p = handles.p;
    udd.col = col;
    udd.nc = length(clusters);
    handles.pop = uicontrol(handles.fi,...
        'style','popupmenu',...
        'position',[20 20 100 40],...
        'string',string,...
        'callback',{@doSelectCluster},...
        'BusyAction','cancel',...
        'Interruptible','off',...
        'visible',visible,'tag',tag);
    set(handles.pop,'units','normalized')
    handles.sli = uicontrol(handles.fi,...
        'style','slider',...
        'position',[50 10 30 20],'max',udd.nc,...
        'sliderstep',[1./(udd.nc+0) 1./(udd.nc+0)],...
        'callback',{@doSwitch2nextCluster},...
        'BusyAction','cancel',...
        'Interruptible','off',...
        'visible',visible,'tag',tag);
    set(handles.sli,'units','normalized')
    udd.pop = handles.pop;
    udd.sli = handles.sli;
    set(handles.pop,'userdata',udd);
    set(handles.sli,'userdata',udd);
end

% Texture thresholding sliders
if  ~isequal(texture,'none') && isequal(clusters,'none')
    if subplotBIN
        subplot(2,1,1)
    end
    udd.tex0 = texture;
    udd.col = col;
    handles.hc = colorbar('peer',ParentAxes);
    set(handles.hc,'visible',visible)
    increment = 0.01;
    % right slider
    handles.s1 = uicontrol(handles.fi,...
        'style','slider',...
        'position',[440 28    20   380],...
        'min',0,'max',length(udd.col),'value',0,...
        'sliderstep',[increment increment],...
        'tooltipstring','texture thresholding control',...
        'callback',{@doThresh},...
        'BusyAction','cancel',...
        'Interruptible','off',...
        'visible',visible,'tag',tag);
    set(handles.s1,'units','normalized')
    udd.s1 = handles.s1;
    % left slider
    handles.s2 = uicontrol(handles.fi,...
        'style','slider',...
        'position',[420 28    20   380],...
        'min',1,'max',length(udd.col),...
        'value',length(udd.col),...
        'sliderstep',[increment increment],...
        'tooltipstring','texture thresholding control',...
        'callback',{@doThresh},...
        'BusyAction','cancel',...
        'Interruptible','off',...
        'visible',visible,'tag',tag);
    set(handles.s2,'units','normalized')
    udd.s2 = handles.s2;
    set(handles.s1,'userdata',udd);
    set(handles.s2,'userdata',udd);
end



set(handles.fi,'visible','on');
drawnow
% if ~addMesh
    try,camlight;end
% end

%cameratoolbar(handles.fi,'setmode','orbit')

out.hfra = getframe(gcf);
out.handles = handles;
out.m = m;

%--------- subfunctions : BUTTONS CALLBACKS ------------%


function doDelMesh(btn,evd)
renderer=get(btn,'userdata');
set(gcf,'renderer',renderer);

function doTransp(btn,evd)
v00=get(btn,'value');
p00=get(btn,'userdata');
set(p00,'facealpha',v00);


function doThresh(btn,evd)
udd00 = get(btn,'userdata');
ind00 = round(get(udd00.s1,'value'));
ind200 = round(get(udd00.s2,'value'));
if(ind200>ind00)
    udd00.col(1:ind00,:)=0.5*ones(ind00,3);
    udd00.col(ind200+1:end,:)=0.5*ones(size(udd00.col(ind200+1:end,:)));
else
    udd00.col(ind200:ind00,:)=0.5*ones(size(udd00.col(ind200:ind00,:)));
end
colormap(udd00.col);
udd00.cax = caxis;


function doSelectCluster(btn,evd)
udd00 = get(btn,'userdata');
ind00=get(gcbo,'value');
set(udd00.sli,'value',ind00-1);
set(udd00.p,'facevertexcdata',udd00.tex(:,ind00));
if ind00 == 1
    colormap(udd00.col);
else
    col00 = colormap(jet);
    col00(1:end/2,:)=0.5*ones(size(col00(1:end/2,:)));
    colormap(col00);
end
udd00.cax = caxis;


function doSwitch2nextCluster(btn,evd)
v00=get(btn,'value')+1;
udd00=get(gcbo,'userdata');
ind00=min([v00 udd00.nc+1]);
set(udd00.pop,'value',ind00);
set(udd00.p,'facevertexcdata',udd00.tex(:,ind00));
if ind00 == 1
    colormap(udd00.col);
else
    col00 = colormap(jet);
    col00(1:end/2,:)=0.5;
    colormap(col00);
end
udd00.cax = caxis;

