function [out] = spm_eeg_displayECD(Pos,Orient,Var,Names,options)
% Plot dipole positions onto the SPM canonical mesh
% FORMAT [out] = spm_eeg_displayDipoles(Pos,Orient,Var,Names,options)
%
% IN (admissible choices):
%   - Pos: a 3xndip matrix containing the positions of the dipoles in
%   the canonical frame of reference
%   - Orient: the same with dipole orientations
%   - Var: the same with position variance
%   - Names: the same with dipole names
%   - options: an optional structure containing
%       .hfig: the handle of the display figure
%       .tag: the tag to be associated with the created UI objects
%       .add: binary variable ({0}, 1: just add dipole in the figure .hfig)
%
% OUT:
%   - out: a structure containing the handles of the object in the figure
%   (including the mesh, the dipoles, the transparency slider, etc...)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_displayECD.m 4054 2010-08-27 19:27:09Z karl $


% checks and defaults
%--------------------------------------------------------------------------
hfig       = [];
ParentAxes = [];
query      = [];
handles    = [];
tag        = '';
try, options; catch, options = [];      end
try, hfig        = options.hfig;        end
try, tag         = options.tag;         end
try, ParentAxes  = options.ParentAxes;  end
try, query       = options.query;       end
try, handles     = options.handles;     end
try
    figure(hfig);
catch
    hfig  = spm_figure('GetWin','Graphics');
    spm_figure('Clear',hfig);
    ParentAxes = axes('parent',hfig);
end
try
    markersize = options.markersize;
catch
    markersize = 20;
end
try
    meshsurf = options.meshsurf;
catch
    meshsurf = fullfile(spm('Dir'),'canonical','cortex_5124.surf.gii');
end

if isscalar(Var), Var = Pos*0 + Var^2;   end
try, Pos{1};    catch, Pos = {Pos};      end
try, Orient{1}; catch, Orient = {Orient};end
try, Var{1};    catch, Var = {Var};      end

ndip = size(Pos{1},2);
if ~exist('Names','var') || isempty(Names)
    for i=1:ndip
        Names{i} = num2str(i);
    end
end


col = ['b','g','r','c','m','y','k','w'];
tmp = ceil(ndip./numel(col));
col = repmat(col,1,tmp);
pa  = get(ParentAxes,'position');

if ndip > 0
    
    if isempty(query)
        opt.hfig = hfig;
        opt.ParentAxes = ParentAxes;
        opt.visible = 'off';
        pos2 = [pa(1),pa(2)+0.25*pa(4),0.03,0.5*pa(4)];
        out  = spm_eeg_render(meshsurf,opt);
        handles.mesh = out.handles.p;
        handles.BUTTONS.transp = out.handles.transp;
        handles.hfig = out.handles.fi;
        handles.ParentAxes  = out.handles.ParentAxes;
        set(handles.mesh,...
            'facealpha',0.1,...
            'visible','on')
        set(handles.BUTTONS.transp,...
            'value',0.1,...
            'position',pos2,...
            'visible','on')
    end
    
    set(ParentAxes,'nextplot','add')
    for j=1:length(Pos)
        for i =1:ndip
            try
                set(handles.hp(j,i),...
                    'xdata',Pos{j}(1,i),...
                    'ydata',Pos{j}(2,i),...
                    'zdata',Pos{j}(3,i));
            catch
                handles.hp(j,i) = plot3(handles.ParentAxes,...
                    Pos{j}(1,i),Pos{j}(2,i),Pos{j}(3,i),...
                    [col(i),'.'],...
                    'markerSize',markersize,...
                    'visible','off');
            end
            try
                no = sqrt(sum(Orient{j}(:,i).^2));
                if no > 0
                    Oi = 1e1.*Orient{j}(:,i)./no;
                else
                    Oi = 1e-5*ones(3,1);
                end
                try
                    set(handles.hq(j,i),...
                        'xdata',Pos{j}(1,i),...
                        'ydata',Pos{j}(2,i),...
                        'zdata',Pos{j}(3,i),...
                        'udata',Oi(1),...
                        'vdata',Oi(2),...
                        'wdata',Oi(3))
                catch
                    handles.hq(j,i) = quiver3(handles.ParentAxes,...
                        Pos{j}(1,i),Pos{j}(2,i),Pos{j}(3,i),...
                        Oi(1),Oi(2),Oi(3),col(i),...
                        'lineWidth',2,'visible','off');
                end
                if isequal(query,'add')
                    set(handles.hq(j,i),...
                        'LineStyle','--',...
                        'lineWidth',1)
                end
            end
            [x,y,z]= ellipsoid(Pos{j}(1,i),Pos{j}(2,i),Pos{j}(3,i),...
                1.*sqrt(Var{j}(1,i)),1.*sqrt(Var{j}(2,i)),1.*sqrt(Var{j}(1,i)),20);
            try
                set(handles.hs(j,i),...
                    'xdata',x,...
                    'ydata',y,...
                    'zdata',z);
            catch
                handles.hs(j,i) = surf(handles.ParentAxes,...
                    x,y,z,...
                    'edgecolor','none',...
                    'facecolor',col(i),...
                    'facealpha',0.2,...
                    'visible','off');
            end
            try
                set(handles.ht(j,i),...
                    'position',Pos{j}(:,i));
            catch
                handles.ht(j,i) = text(...
                    Pos{j}(1,i),Pos{j}(2,i),Pos{j}(3,i),...
                    Names{i},...
                    'Parent',handles.ParentAxes,...
                    'visible','off');
            end
        end
    end
    
    if length(Pos) > 1
        
        try, set(handles.hp(end,:),'visible','on'); end
        try, set(handles.hq(end,:),'visible','on'); end
        try, set(handles.hs(end,:),'visible','on'); end
        try, set(handles.ht(end,:),'visible','on'); end
        
        handles.uic(1) = uicontrol(handles.fi,...
            'units','normalized',...
            'position',[0.45,0.5,0.2,0.03],...
            'style','radio','string','Show priors',...
            'callback',@doChange1,...
            'BackgroundColor',[1 1 1],...
            'tooltipstring','Display prior locations',...
            'userdata',handles,'value',0,...
            'BusyAction','cancel',...
            'Interruptible','off',...
            'tag','plotEEG');

        handles.uic(2) = uicontrol(handles.fi,...
            'units','normalized',...
            'position',[0.45,0.53,0.2,0.03],...
            'style','radio','string','Show posteriors',...
            'callback',@doChange2,...
            'BackgroundColor',[1 1 1],...
            'tooltipstring','Display posterior locations',...
            'userdata',handles,'value',1,...
            'BusyAction','cancel',...
            'Interruptible','off',...
            'tag','plotEEG');
        
    else
        
        try, set(handles.hp(1,:),'visible','on'); end
        try, set(handles.hq(1,:),'visible','on'); end
        try, set(handles.hs(1,:),'visible','on'); end
        try, set(handles.ht(1,:),'visible','on'); end
        
    end
    
end

try
    clear out
    out.handles = handles;
catch
    out = [];    
end

%==========================================================================
function doChange1(i1,i2)
val = get(i1,'value');
handles = get(i1,'userdata');
if ~val
    try, set(handles.hp(1,:),'visible','off'); end
    try, set(handles.hq(1,:),'visible','off'); end
    try, set(handles.hs(1,:),'visible','off'); end
    try, set(handles.ht(1,:),'visible','off'); end
else
    try, set(handles.hp(1,:),'visible','on');  end
    try, set(handles.hq(1,:),'visible','on');  end
    try, set(handles.hs(1,:),'visible','on');  end
    try, set(handles.ht(1,:),'visible','on');  end
end


%==========================================================================
function doChange2(i1,i2)
val = get(i1,'value');
handles = get(i1,'userdata');
if ~val
    try, set(handles.hp(2,:),'visible','off'); end
    try, set(handles.hq(2,:),'visible','off'); end
    try, set(handles.hs(2,:),'visible','off'); end
    try, set(handles.ht(2,:),'visible','off'); end
else
    try, set(handles.hp(2,:),'visible','on');  end
    try, set(handles.hq(2,:),'visible','on');  end
    try, set(handles.hs(2,:),'visible','on');  end
    try, set(handles.ht(2,:),'visible','on');  end
end
