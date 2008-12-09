function  [out] = spm_eeg_render(m,options)
%
% [out] = spm_eeg_render(m,options)
%
% This function is a visualization routine, mainly for texture and
% clustering on the cortical surface.
% IN:
% - m = MATLAB mesh (containing the fields .faces et .vertices)
% - options = structure variable:
%       .texture = texture to be projected onto the mesh
%       .clusters = cortical parcelling (cell variable containing the
%       vertex indices of each cluster) 
%       .clustersName = name of the clusters
%       .figname = name to be given to the window
%       .ParentAxes = handle of the axes within which the mesh should be
%       displayed
%       .hfig = handle of existing figure. If this option is provided, then
%       visu_maillage_surf adds the (textured) mesh to the figure hfig, and
%       a control for its transparancy.
%       .subplot.bin = flag for figure dividing
%       .subplot.y = 2D matrix to be plotted in a subfigure (only if the
%       associated options.subplot.bin = 1)
% OUT:
%   - out: a structure containing the fields:
%       .hfra: frame structure for movie building
%       .handles: a structure containing the handles of the created
%       uicontrols and mesh objects.
% NB: The texture and the clusters can not be visualized at the same time.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_render.m 2540 2008-12-09 17:14:25Z jean $




%----------------------------------------------------------------------%
%------------- Common features for any visualization ------------------%
%----------------------------------------------------------------------%

% Default options
handles.fi = figure('visible','off',...
    'NumberTitle','Off','Name','Mesh visualization',...
    'tag','visu_maillage_surf');
ns = 0;
texture = 'none';
clusters = 'none';
figname = 'none';
subplotBIN = 0;
addMesh = 0;
tag = '';
visible = 'on';
ParentAxes = gca;

try,options;catch,options = [];end

% Now get options
if ~isempty(options)
    
    % get texture if provided
    if isfield(options,'texture')
        texture = options.texture;
    end
    
    if isfield(options,'ParentAxes') && ~isempty(options.ParentAxes)
        ParentAxes = options.ParentAxes;
    end
    
    if isfield(options,'tag')
        tag = options.tag;
    end

    if isfield(options,'visible')
        visible = options.visible;
    end
    
    % get custers if provided
    if isfield(options,'clusters')
        figname = 'Parcelling vizualization';
        clusters = options.clusters;
        IND = zeros(1,length(m.vertices));
        K = length(clusters);
        for k = 1:K
            IND(clusters{k}) = k+1./K;
        end
        texture = IND';
    end

    % get figname if provided
    if isfield(options,'figname')
        figname = options.figname;
        set(handles.fi,'NumberTitle','Off','Name',figname);
    end
    
    % get figure handle for possible multiplot
    if isfield(options,'hfig')
        try
            figure(options.hfig)
            ParentAxes = gca;
            hold on
            close(handles.fi);
            handles.fi = options.hfig;
            addMesh = 1;
        end
        try % get number of transparency sliders in current figure...
            hh=get(handles.fi,'children');
            ns=length(findobj(hh,'userdata','tag_UIC_transparency'));
        catch
            ns=1;
        end
    end
    
    % Subdivide figure?
    if isfield(options,'subplot') && options.subplot.bin
        subplotBIN = 1;
        subplot(2,1,1)
    end
    
else
   
    
end

oldRenderer = get(handles.fi,'renderer');
set(handles.fi,'renderer','OpenGL');


% Plot mesh and texture/clusters
if isequal(texture,'none') == 1
    figure(handles.fi)
    handles.p = patch(m);
    set(handles.p, 'facecolor', [.5 .5 .5], 'EdgeColor', 'none',...
        'parent',ParentAxes,...
        'userdata',oldRenderer,...
        'visible',visible,'tag',tag);
else
    texture = texture(:);
    figure(handles.fi)
    if isequal(length(texture),length(m.vertices))
        handles.p = patch(m,'facevertexcdata',texture,...
            'facecolor','interp','EdgeColor', 'none',...
            'parent',ParentAxes,...
            'userdata',oldRenderer,...
            'visible',visible,'tag',tag);
        colormap(jet(256));
        udd.tex = texture;
        udd.cax = caxis;
    else
        texture = 'none';
        disp('Warning: size of texture does not match number of vertices!')
        handles.p = patch(m,'facecolor', [.5 .5 .5], 'EdgeColor', 'none',...
            'parent',ParentAxes,...
            'userdata',oldRenderer,...
            'visible',visible,'tag',tag);
    end
end

set(handles.p,'deleteFcn',@doDelMesh)


figure(handles.fi)  % avoid user confusion with the open windows
daspect([1 1 1]);
view(142.5,30);
axis tight;
axis off
lighting gouraud;
camva('auto');
view(25,45);


% plot second subplot if provided
if subplotBIN
    figure(handles.fi)
    subplot(2,1,1),zoom(2);
    subplot(2,1,2),plot(options.subplot.y')
    set(gca,'ygrid','on')
    set(gca,'ylim',[-0.2,1.2])
    set(gca,'box','off')
end


% build internal userdata structure
udd.p = handles.p;

% xlabel('x')
% ylabel('y')
% zlabel('z')








%----------------------------------------------------------------------%
%---------------------- GUI tools and buttons -------------------------%
%----------------------------------------------------------------------%


% Transparancy sliders
pos = [20 100 20 245];
pos(1) = pos(1) + ns.*25;
handles.transp = uicontrol(...
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
handles.tag = uicontrol(...
    'style','text',...
    'visible','off',...
    'tag',tag,...
    'userdata','tag_UIC_transparency');
set(handles.transp,'units','normalized')
udd.transp = handles.transp;
% 
% % Cortex inflation togglebutton
% posh = [200 20 100 40];
% hu=uicontrol('style','togglebutton','position',posh);
% set(hu,'string','inflate cortex');
% ud2.m = m;
% ud2.p = p;
% ud2.texture = texture;
% set(hu,'callback','vms_inflate')
% set(hu,'userdata',ud2);
% set(hu,'enable','off')  % still to be debugged
% set(hu,'units','normalized')


% 
% % SpitBrain checkboxes
% pos = [20 360 80 20];
% sb1 = uicontrol('style','checkbox','position',pos,'callback','vms_splitBrain');
% set(sb1,'string','left hemi','value',1)
% set(sb1,'units','normalized')
% pos = [20 390 80 20];
% sb2 = uicontrol('style','checkbox','position',pos,'callback','vms_splitBrain');
% set(sb2,'string','right hemi','value',1)
% set(sb2,'units','normalized')
% udd.sb1 = sb1;
% udd.sb2 = sb2;
% udd.m = m;
% udd.split = 0;
% udd.fi = fi;



% Clustering buttons and popup menu
if ~isequal(clusters,'none')
    if subplotBIN
        subplot(2,1,1)
    end
    %             set(p,'FaceColor','flat');
    col=colormap(lines);
    nc = floor(256./K);
    col = [repmat([0.8157    0.6666    0.5762],nc/2,1);kron(col(1:K,:),ones(nc,1))];
    if K > 1
        col(end-nc/2:end,:) = [];
    end
    colormap(col);
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
    handles.pop = uicontrol('style','popupmenu',...
        'position',[20 20 100 40],...
        'string',string,...
        'callback',{@doSelectCluster},...
        'BusyAction','cancel',...
        'Interruptible','off',...
        'visible',visible,'tag',tag);
    set(handles.pop,'units','normalized')
    handles.sli = uicontrol('style','slider',...
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
    tmp = colormap;
    udd.col = tmp;
    handles.hc = colorbar;
    set(handles.hc,'visible',visible);
    increment = 0.01;
    % right slider
    handles.s1 = uicontrol('style','slider',...
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
    handles.s2 = uicontrol('style','slider',...
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

% 
% set(handles.sb1,'userdata',udd);
% set(handles.sb2,'userdata',udd);


set(handles.fi,'visible','on');
drawnow
% if ~addMesh
    camlight
% end

cameratoolbar(handles.fi,'setmode','orbit')


gco;

out.hfra = getframe(gcf);
out.handles = handles;


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
% set(udd00.sb1,'userdata',udd00);
% set(udd00.sb2,'userdata',udd00);


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
% set(udd00.sb1,'userdata',udd00);
% set(udd00.sb2,'userdata',udd00);


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
    col00(1:end/2,:)=0.5;%*ones(size(col00(1:end/2,:)));
    colormap(col00);
end
udd00.cax = caxis;
% set(udd00.sb1,'userdata',udd00);
% set(udd00.sb2,'userdata',udd00);

