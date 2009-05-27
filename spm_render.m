function spm_render(dat,brt,rendfile)
% Render blobs on surface of a 'standard' brain
% FORMAT spm_render(dat,brt,rendfile)
%
% dat      - a struct array of length 1 to 3
%            each element is a structure containing:
%            - XYZ - the x, y & z coordinates of the transformed SPM{.}
%                    values in units of voxels.
%            - t   - the SPM{.} values.
%            - mat - affine matrix mapping from XYZ voxels to MNI.
%            - dim - dimensions of volume from which XYZ is drawn.
% brt      - brightness control:
%            If NaN, then displays using the old style with hot
%            metal for the blobs, and grey for the brain.
%            Otherwise, it is used as a ``gamma correction'' to
%            optionally brighten the blobs up a little.
% rendfile - the file containing the images to render on to (see also
%            spm_surf.m) or a surface mesh file.
%
% Without arguments, spm_render acts as its own UI.
%__________________________________________________________________________
% 
% spm_render prompts for details of up to three SPM{.}s that are then
% displayed superimposed on the surface of a 'standard' brain.
%
% The first is shown in red, then green then blue.
%
% The blobs which are displayed are the integral of all transformed t
% values, exponentially decayed according to their depth. Voxels that
% are 10mm behind the surface have half the intensity of ones at the
% surface.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_render.m 3152 2009-05-27 10:54:49Z guillaume $

SVNrev = '$Rev: 3152 $';

global prevrend
if ~isstruct(prevrend)
    prevrend = struct('rendfile','',...
                      'brt',[],...
                      'col',[]);
end

%-Parse arguments, get data if not passed as parameters
%==========================================================================
if nargin < 1
    spm('FnBanner',mfilename,SVNrev);
    spm('FigName','Results: render');

    num   = spm_input('Number of sets',1,'1 set|2 sets|3 sets',[1 2 3],1);

    for i = 1:num
        [SPM,xSPM] = spm_getSPM;
        dat(i)    = struct( 'XYZ',  xSPM.XYZ,...
                    't',    xSPM.Z',...
                    'mat',  xSPM.M,...
                    'dim',  xSPM.DIM);
    end
    showbar = 1;
else
    num     = length(dat);
    showbar = 0;
end

%-Get surface
%--------------------------------------------------------------------------
if nargin < 3 || isempty(prevrend.rendfile)
    [rendfile, sts] = spm_select(1,'mesh','Render file'); % .mat or .gii file
    if ~sts, return; end
end
prevrend.rendfile = rendfile;

%-Get brightness & colours
%--------------------------------------------------------------------------
if nargin < 2  || isempty(prevrend.brt)
    brt = 1;
    if num==1
        brt = spm_input('Style',1,'new|old',[1 NaN], 1);
    end

    if isfinite(brt)
        brt = spm_input('Brighten blobs',1,'none|slightly|more|lots',[1 0.75 0.5 0.25], 1);
        col = eye(3);
        % ask for custom colours & get rgb values
        %------------------------------------------------------------------
        if spm_input('Which colours?','!+1','b',{'RGB','Custom'},[0 1],1)
            for k = 1:num,
                col(k,:) = uisetcolor(col(k,:),sprintf('Color of blob set %d',k));
            end
        end
    else
        col = [];
    end
elseif isfinite(brt) && isempty(prevrend.col)
    col = eye(3);
elseif isfinite(brt)  % don't need to check prevrend.col again
    col = prevrend.col;
else
    col = [];
end
prevrend.brt = brt;
prevrend.col = col;

%-Perform the rendering
%==========================================================================
[p,f,e] = fileparts(rendfile);
loadgifti = false;
if strcmpi(e,'.mat')
    load(rendfile);
    if ~exist('rend','var') && ~exist('Matrixes','var')
        loadgifti = true;
    end
end
if ~strcmpi(e,'.mat') || loadgifti
    try
        rend = export(gifti(rendfile),'patch');
    catch
        error('\nCannot read  render file "%s".\n', rendfile);
    end
    surf_rend(dat,rend,col);
    return
end

spm('Pointer','Watch');

if ~exist('rend','var') % Assume old format...
    rend = cell(size(Matrixes,1),1);
    for i=1:size(Matrixes,1),
        rend{i}=struct('M',eval(Matrixes(i,:)),...
            'ren',eval(Rens(i,:)),...
            'dep',eval(Depths(i,:)));
        rend{i}.ren = rend{i}.ren/max(max(rend{i}.ren));
    end
end

if showbar, spm_progress_bar('Init', size(dat,1)*length(rend),...
            'Formatting Renderings', 'Number completed'); end
for i=1:length(rend),
    rend{i}.max=0;
    rend{i}.data = cell(size(dat,1),1);
    if issparse(rend{i}.ren),
        % Assume that images have been DCT compressed
        % - the SPM99 distribution was originally too big.
        d = size(rend{i}.ren);
        B1 = spm_dctmtx(d(1),d(1));
        B2 = spm_dctmtx(d(2),d(2));
        rend{i}.ren = B1*rend{i}.ren*B2';
        % the depths did not compress so well with
        % a straight DCT - therefore it was modified slightly
        rend{i}.dep = exp(B1*rend{i}.dep*B2')-1;
    end
    rend{i}.ren(rend{i}.ren>=1) = 1;
    rend{i}.ren(rend{i}.ren<=0) = 0;
    if showbar, spm_progress_bar('Set', i); end
end
if showbar, spm_progress_bar('Clear'); end

if showbar, spm_progress_bar('Init', length(dat)*length(rend),...
            'Making pictures', 'Number completed'); end

mx = zeros(length(rend),1)+eps;
mn = zeros(length(rend),1);

for j=1:length(dat),
    XYZ = dat(j).XYZ;
    t   = dat(j).t;
    dim = dat(j).dim;
    mat = dat(j).mat;

    for i=1:length(rend),

        % transform from Talairach space to space of the rendered image
        %------------------------------------------------------------------
        M1  = rend{i}.M*mat;
        zm  = sum(M1(1:2,1:3).^2,2).^(-1/2);
        M2  = diag([zm' 1 1]);
        M  = M2*M1;
        cor = [1 1 1 ; dim(1) 1 1 ; 1 dim(2) 1; dim(1) dim(2) 1 ;
               1 1 dim(3) ; dim(1) 1 dim(3) ; 1 dim(2) dim(3); dim(1) dim(2) dim(3)]';
        tcor= M(1:3,1:3)*cor + M(1:3,4)*ones(1,8);
        off = min(tcor(1:2,:)');
        M2  = spm_matrix(-off+1)*M2;
        M  = M2*M1;
        xyz = (M(1:3,1:3)*XYZ + M(1:3,4)*ones(1,size(XYZ,2)));
        d2  = ceil(max(xyz(1:2,:)'));

        % Calculate 'depth' of values
        %------------------------------------------------------------------
        dep = spm_slice_vol(rend{i}.dep,spm_matrix([0 0 1])*inv(M2),d2,1);
        z1  = dep(round(xyz(1,:))+round(xyz(2,:)-1)*size(dep,1));

        if ~isfinite(brt), msk = find(xyz(3,:) < (z1+20) & xyz(3,:) > (z1-5));
        else,      msk = find(xyz(3,:) < (z1+60) & xyz(3,:) > (z1-5)); end

        if ~isempty(msk),

            % Generate an image of the integral of the blob values.
            %--------------------------------------------------------------
            xyz = xyz(:,msk);
            if ~isfinite(brt), t0  = t(msk);
            else,   dst = xyz(3,:) - z1(msk);
                dst = max(dst,0);
                t0  = t(msk).*exp((log(0.5)/10)*dst)';
            end
            X0  = full(sparse(round(xyz(1,:)), round(xyz(2,:)), t0, d2(1), d2(2)));
            hld = 1; if ~isfinite(brt), hld = 0; end
            X   = spm_slice_vol(X0,spm_matrix([0 0 1])*M2,size(rend{i}.dep),hld);
            msk = find(X<0);
            X(msk) = 0;
        else
            X = zeros(size(rend{i}.dep));
        end

        % Brighten the blobs
        %------------------------------------------------------------------
        if isfinite(brt), X = X.^brt; end

        mx(j) = max([mx(j) max(max(X))]);
        mn(j) = min([mn(j) min(min(X))]);

        rend{i}.data{j} = X;

        if showbar, spm_progress_bar('Set', i+(j-1)*length(rend)); end
    end
end

mxmx = max(mx);
mnmn = min(mn);

if showbar, spm_progress_bar('Clear'); end
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);

nrow = ceil(length(rend)/2);
if showbar, hght = 0.95; else, hght = 0.5; end
% subplot('Position',[0, 0, 1, hght]);
ax=axes('Parent',Fgraph,'units','normalized','Position',[0, 0, 1, hght],'Visible','off');
image(0,'Parent',ax);
set(ax,'YTick',[],'XTick',[]);

if ~isfinite(brt),
    % Old style split colourmap display.
    %----------------------------------------------------------------------
    load Split;
    colormap(split);
    for i=1:length(rend),
        ren = rend{i}.ren;
        X   = (rend{i}.data{1}-mnmn)/(mxmx-mnmn);
        msk = find(X);
        ren(msk) = X(msk)+(1+1.51/64);
        ax=axes('Parent',Fgraph,'units','normalized',...
            'Position',[rem(i-1,2)*0.5, floor((i-1)/2)*hght/nrow, 0.5, hght/nrow],...
            'Visible','off');
        image(ren*64,'Parent',ax);
        set(ax,'DataAspectRatio',[1 1 1], ...
            'PlotBoxAspectRatioMode','auto',...
            'YTick',[],'XTick',[],'XDir','normal','YDir','normal');
    end
else
    % Combine the brain surface renderings with the blobs, and display using
    % 24 bit colour.
    %----------------------------------------------------------------------
    for i=1:length(rend),
        ren = rend{i}.ren;
        X = cell(3,1);
        for j=1:length(rend{i}.data),
            X{j} = rend{i}.data{j}/(mxmx-mnmn)-mnmn;
        end
        for j=(length(rend{i}.data)+1):3
            X{j}=zeros(size(X{1}));
        end

        rgb = zeros([size(ren) 3]);
        tmp = ren.*max(1-X{1}-X{2}-X{3},0);
        for k = 1:3
            rgb(:,:,k) = tmp + X{1}*col(1,k) + X{2}*col(2,k) +X{3}*col(3,k);
        end
        rgb(rgb>1) = 1;         
        
        ax=axes('Parent',Fgraph,'units','normalized',...
            'Position',[rem(i-1,2)*0.5, floor((i-1)/2)*hght/nrow, 0.5, hght/nrow],...
            'nextplot','add', ...
            'Visible','off');
        image(rgb,'Parent',ax);
        set(ax,'DataAspectRatio',[1 1 1], ...
            'PlotBoxAspectRatioMode','auto',...
            'YTick',[],'XTick',[],...
            'XDir','normal','YDir','normal');
    end
end

spm('Pointer','Arrow');

%==========================================================================
% function surf_rend(dat,rend,col)
%==========================================================================
function surf_rend(dat,rend,col)

Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);
rdr = get(Fgraph,'Renderer');
set(Fgraph,'Renderer','OpenGL');

ax = axes(...
    'Parent',Fgraph,...
    'units','normalized',...
    'Position',[0.025, 0.025, 0.95, 0.45],...
    'Visible','off');

Vo = spm_write_filtered(dat.t,dat.XYZ,dat.dim,dat.mat,'',tempname);
v = spm_get_data(Vo,...
    double(inv(Vo.mat)*[rend.vertices';ones(1,size(rend.vertices,1))]));
spm_unlink(Vo.fname); spm_unlink([spm_str_manip(Vo.fname,'s') '.hdr']);
C = spm_mesh_curvature(rend) > 0;
C = 0.5 * repmat(C,1,3) + 0.3 * repmat(~C,1,3);
col = hot(256);
col(1,:) = 0.5;
if ~any(v)
    cdat = 0.5*ones(length(v),3);
else
    cdat = squeeze(ind2rgb(floor(v(:)/max(v(:))*size(col,1)),col));
end
cdat = repmat(v==0,3,1)' .* C + repmat(v~=0,3,1)' .* cdat;

hp = patch(rend, 'Parent',ax,...
    'FaceVertexCData',cdat, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none',...
    'FaceLighting', 'phong',...
    'SpecularStrength' ,0.7, 'AmbientStrength', 0.1,...
    'DiffuseStrength', 0.7, 'SpecularExponent', 10);

set(Fgraph,'CurrentAxes',ax);
view(ax,[-90 0]);
axis(ax,'image');

l = camlight; set(l,'Parent',ax);
material(Fgraph,'dull');
setappdata(ax,'camlight',l);

%spm_mesh_inflate(hp,Inf,1);
%view(ax,[-90 0]);axis(ax,'image');

r = rotate3d(ax);
set(r,'enable','off');
cmenu = uicontextmenu;
c1 = uimenu(cmenu, 'Label', 'Inflate', 'Interruptible','off', 'Callback', @myinflate);
setappdata(c1,'patch',hp);
setappdata(c1,'axis',ax);
c12 = uimenu(cmenu, 'Label', 'Connected Components', 'Visible', 'off', 'Interruptible','off', 'Callback', @mycclabel);
C = spm_mesh_label(hp);
setappdata(c12,'patch',hp);
setappdata(c12,'cclabel',C);
c2 = uimenu(cmenu, 'Label', 'Rotate', 'Checked','on','Separator','on','Callback', @myswitchrotate);
setappdata(c2,'rotate3d',r);
c3 = uimenu(cmenu, 'Label', 'Transparency');
setappdata(c3,'patch',hp);
uimenu(c3,'Label','0%',  'Checked','on',  'Callback',@mytransparency);
uimenu(c3,'Label','20%', 'Checked','off', 'Callback',@mytransparency);
uimenu(c3,'Label','40%', 'Checked','off', 'Callback',@mytransparency);
uimenu(c3,'Label','60%', 'Checked','off', 'Callback',@mytransparency);
uimenu(c3,'Label','80%', 'Checked','off', 'Callback',@mytransparency);
c4 = uimenu(cmenu, 'Label', 'Save As...','Separator','on','Callback', @mysave);
setappdata(c4,'patch',hp);
try, set(r,'uicontextmenu',cmenu); end
try, set(hp,'uicontextmenu',cmenu); end
set(r,'enable','on');
set(r,'ActionPostCallback',@mypostcallback);
set(hp,'DeleteFcn',{@mydeletefcn,Fgraph,rdr});

try % meaningless when called outside spm_results_ui
    hReg = spm_XYZreg('FindReg',spm_figure('GetWin','Interactive'));
    xyz  = spm_XYZreg('GetCoords',hReg);
    hs   = mydispcursor('Create',ax,dat.mat,xyz);
    spm_XYZreg('Add2Reg',hReg,hs,@mydispcursor);
end
    
%==========================================================================
function myinflate(obj,evd)
spm_mesh_inflate(getappdata(obj,'patch'),Inf,1);
axis(getappdata(obj,'axis'),'image');

%==========================================================================
function mycclabel(obj,evd)
C = getappdata(obj,'cclabel');
F = get(getappdata(obj,'patch'),'Faces');
V = get(getappdata(obj,'patch'),'Vertices');
V = zeros(size(V,1),1);
V(reshape(F(C==1,:),[],1)) = 1;
set(getappdata(obj,'patch'),'FaceVertexAlphaData',V);
set(getappdata(obj,'patch'),'FaceAlpha','interp');

%==========================================================================
function myswitchrotate(obj,evd)
if strcmpi(get(getappdata(obj,'rotate3d'),'enable'),'on')
    set(getappdata(obj,'rotate3d'),'enable','off');
    set(obj,'Checked','off');
else
    set(getappdata(obj,'rotate3d'),'enable','on');
    set(obj,'Checked','on');
end

%==========================================================================
function mytransparency(obj,evd)
t = 1 - sscanf(get(obj,'Label'),'%d%%') / 100;
set(getappdata(get(obj,'parent'),'patch'),'FaceAlpha',t);
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');

%==========================================================================
function mysave(obj,evd)
g = gifti;
g.vertices = get(getappdata(obj,'patch'),'Vertices');
g.faces = get(getappdata(obj,'patch'),'Faces');
g.cdata = get(getappdata(obj,'patch'),'FaceVertexCData');
[filename, pathname] = uiputfile({'*.gii' 'GIfTI files (*.gii)'}, 'Save as');
if ~isequal(filename,0) && ~isequal(pathname,0)
    [pth,nam,ext] = fileparts(filename);
    if ~strcmpi(ext,'.gii'), filename = [filename '.gii']; end
	save(g,fullfile(pathname, filename));
end

%==========================================================================
function mypostcallback(obj,evd)
try, camlight(getappdata(evd.Axes,'camlight')); end

%==========================================================================
function mydeletefcn(obj,evd,varargin)
try, rotate3d(get(obj,'parent'),'off'); end
set(varargin{1},'Renderer',varargin{2});

%==========================================================================
function varargout = mydispcursor(varargin)

switch lower(varargin{1})
    %======================================================================
    case 'create'
    %======================================================================
    % hMe = mydispcursor('Create',ax,M,xyz)
    ax  = varargin{2};
    M   = varargin{3};
    xyz = varargin{4};
    
    [X,Y,Z] = sphere;
    vx = sqrt(sum(M(1:3,1:3).^2));
    X = X*vx(1) + xyz(1);
    Y = Y*vx(2) + xyz(2);
    Z = Z*vx(3) + xyz(3);
    hold(ax,'on');
    hs = surf(X,Y,Z,'parent',ax,...
        'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting', 'phong');
    set(hs,'UserData',xyz);
    
    varargout = {hs};
    
    %=======================================================================
    case 'setcoords'    % Set co-ordinates
    %=======================================================================
    % [xyz,d] = mydispcursor('SetCoords',xyz,hMe,hC)
    hMe  = varargin{3};
    pxyz = get(hMe,'UserData');
    xyz  = varargin{2};
    
    set(hMe,'XData',get(hMe,'XData') - pxyz(1) + xyz(1));
    set(hMe,'YData',get(hMe,'YData') - pxyz(2) + xyz(2));
    set(hMe,'ZData',get(hMe,'ZData') - pxyz(3) + xyz(3));
    set(hMe,'UserData',xyz);
    
    varargout = {xyz,[]};
    
    %=======================================================================
    otherwise
    %=======================================================================
    error('Unknown action string')

end
