function [f] = spm_opm_plotScalpData(S)
% Display M/EEG interpolated sensor data on a scalp image
% FORMAT D = spm_opm_amm(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.T             - time point to initalise to                    - Default: first sample  
%   S.display       - string to deermine what is plotted   -Default: 'RADIAL'
% OUTPUT:
%   f          - the handle of the figure which displays the interpolated
%                data
%__________________________________________________________________________
%
% This function creates a figure whose purpose is to display an
% interpolation of the sensor data on the scalp (as an image).
%__________________________________________________________________________

% Jean Daunizeau
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
if ~isfield(S, 'display'),        S.display = 'RADIAL'; end
if ~isfield(S, 'noButtons'),      S.noButtons = 0; end


%- Select radial
%--------------------------------------------------------------------------
s = sensors(S.D,'MEG');
pos = s.coilpos;
ori = s.coilori;


[C,~,~]= fitSphere(pos(:,1),pos(:,2),pos(:,3));
cpos = bsxfun(@minus,pos,C);
norm = sqrt(sum(cpos.^2,2));
u = cpos./repmat(norm,1,3);
radial = abs(sum(ori.*u,2))>.8;
notrad = pos(~radial,:);
pos = pos(radial,:);
ChanLabel  = char(s.label(radial));
Z=S.D(indchannel(S.D,cellstr(ChanLabel)),:,:);
notLabel  = char(s.label(~radial));
notZ=S.D(indchannel(S.D,cellstr(notLabel)),:,:);

%- Select norm
%--------------------------------------------------------------------------
if (strmatch(S.display,'NORM'))
  for i = 1:size(pos,1)
    otherChan =sum(abs(bsxfun(@minus,notrad,pos(i,:))),2)<.001;
    Z(i,:)=sqrt(Z(i,:).^2+sum(notZ(otherChan,:).^2,1));
  end
end

%- infer axis
%--------------------------------------------------------------------------

%- regex
%--------------------------------------------------------------------------




ParentAxes = [];
f          = [];
clim       = [min(Z(:))-( max(Z(:))-min(Z(:)) )/63 , max(Z(:))];
figName    = 'Image Scalp data';
noButtons  = S.noButtons;
in     = [];
in.cbar = 1;
in.plotpos = 1;
try 
c2d = S.D.coor2D;
ps= c2d(:,radial);
catch
% get 2D positions from 3D positions
[ps] = get2Dfrom3D(pos);
ps   = ps';
end

% exclude channels ?
goodChannels = find(~isnan(ps(1,:)));
ps          = ps(:,goodChannels);
Z            = Z(goodChannels,:);
ChanLabel    = ChanLabel(goodChannels, :);
cZ         = Z;
cpos       = ps;
cChanLabel = ChanLabel;

xmin    = min(cpos(1,:));
xmax    = max(cpos(1,:));
dx      = (xmax-xmin)./250;
ymin    = min(cpos(2,:));
ymax    = max(cpos(2,:));
dy      = (ymax-ymin)./250;
x       = xmin:dx:xmax;
y       = ymin:dy:ymax;
[XI,YI] = meshgrid(x,y);

if(size(cZ,2)>1 && isfield(S,'T'))
[~,ind]=min(abs(S.D.time-S.T));
  cZ = Z(:,ind);
  T = S.D.time(ind);
else 
  cZ = Z(:,1);
  T = S.D.time(1);
end

ZI      = griddata(cpos(1,:)',cpos(2,:)',full(double(cZ')),XI,YI);

try
    figure(f)
catch
    f   = figure(...
        'name',      figName,...
        'color',     [1 1 1],...
        'deleteFcn', @dFcn);
    ParentAxes = axes('parent',f);    
end

COLOR = get(f,'color');
d.hi = image(flipud(ZI),...
    'CDataMapping','scaled',...
    'Parent',ParentAxes);
set(ParentAxes,'nextPlot','add',...
    'tag','spm_eeg_plotScalpData')
try
    if length(unique(ZI)) ~= 1
        [C,d.hc] = contour(ParentAxes,flipud(ZI),...
            'linecolor',0.5.*ones(3,1));
    end
end
caxis(ParentAxes,clim);
col = jet;
col(1,:) = COLOR;
colormap(ParentAxes,col)

if in.cbar
    d.cbar = colorbar('peer',ParentAxes);
end

axis(ParentAxes,'off')
axis(ParentAxes,'equal')
axis(ParentAxes,'tight')

fpos = cpos;
fpos(1,:) = fpos(1,:) - xmin;
fpos(2,:) = fpos(2,:) - ymin;
fpos(1,:) = fpos(1,:)./(dx);
fpos(2,:) = fpos(2,:)./(dy);
fpos(2,:) = 250-fpos(2,:);  % for display purposes (flipud imagesc)

figure(f);
if in.plotpos
    d.hp = plot(ParentAxes,...
        fpos(1,:),fpos(2,:),...
        'ko');
end

d.ht = text(fpos(1,:),fpos(2,:),cChanLabel,...
    'Parent',ParentAxes,...
    'visible','off');
axis(ParentAxes,'image')

d.interp.XI = XI;
d.interp.YI = YI;
d.interp.pos = cpos;
d.f = f;
d.pos = fpos;
d.goodChannels = goodChannels;
d.ChanLabel = cChanLabel;
d.origChanLabel = ChanLabel;
d.origpos = pos;
d.ParentAxes = ParentAxes;
d.in = in;


if ~noButtons
    d.hsp = uicontrol(f,...
        'style','pushbutton',...
        'callback',{@dosp},...
        'BusyAction','cancel',...
        'Interruptible','off',...
        'position',[10    50    80    20],...
        'string','channel pos');
    d.hsn = uicontrol(f,...
        'style','pushbutton',...
        'callback',{@dosn},...
        'BusyAction','cancel',...
        'Interruptible','off',...
        'position',[10    80    80    20],...
        'string','channel names');
end
if ~isempty(in) && isfield(in,'handles')  
    nT = length(in.gridTime);
    d.hti = uicontrol(f,...
        'style','text',...
        'BackgroundColor',COLOR,...
        'string',[num2str(in.gridTime(in.x)),' (',in.unit,')'],...
        'position',[10    10    120    20]);
    d.hts = uicontrol(f,...
        'style','slider',...
        'Position',[130 10 250 20],...
        'min',1,'max',nT,...
        'value',in.x,'sliderstep',[1./(nT-1) 1./(nT-1)],...
        'callback',{@doChangeTime},...
        'BusyAction','cancel',...
        'Interruptible','off');
    set(d.hti,'userdata',d);
    set(d.hts,'userdata',d);
end

if size(Z,2)>1 && S.noButtons~=1
 nT = size(Z,2);
 in.x = 1;
 d.hti = uicontrol(f,...
        'style','text',...
        'BackgroundColor',COLOR,...
        'string',[num2str(T),' (','s',')'],...
        'position',[10    10    120    20]);
     d.hts = uicontrol(f,...
        'style','slider',...
        'Position',[130 10 250 20],...
        'min',1,'max',nT,...
        'value',in.x,'sliderstep',[4./S.D.fsample 1./S.D.fsample],...
        'callback',{@doChangeTime},...
        'BusyAction','cancel',...
        'Interruptible','off');
     ze = zeros([size(Z),1]);
     ze(:,:,1)=Z;   
    set(f,'userdata',ze);
    d.in.handles.hfig=f;
    d.in.ind = 1:size(Z,1);
    in.gridTime=[];
    if(isempty(in.gridTime))
        d.in.gridTime = S.D.time();
    end
    d.in.type = 'MEG';
    d.in.unit = 's';
    
    set(d.hti,'userdata',d);
    set(d.hts,'userdata',d);
   
end 

if ~noButtons
    set(d.hsp,'userdata',d);
    set(d.hsn,'userdata',d);
end
set(d.ParentAxes,'userdata',d);
if (~strmatch(S.display,'NORM'))
  caxis([-max(abs(Z(:))) max(abs(Z(:)))])
end

end

%==========================================================================
% dFcn
%==========================================================================
function dFcn(btn,evd)
hf = findobj('tag','Graphics');
D = get(hf,'userdata');
try delete(D.PSD.handles.hli); end
end

%==========================================================================
% dosp
%==========================================================================
function dosp(btn,evd)
d = get(btn,'userdata');
switch get(d.hp,'visible')
    case 'on'
        set(d.hp,'visible','off');
    case 'off'
        set(d.hp,'visible','on');
end
end

%==========================================================================
% dosn
%==========================================================================
function dosn(btn,evd)
d = get(btn,'userdata');
switch get(d.ht(1),'visible')
    case 'on'
        set(d.ht,'visible','off');
    case 'off'
        set(d.ht,'visible','on');
end
end

%==========================================================================
% doChangeTime
%==========================================================================
function doChangeTime(btn,evd)
d = get(btn,'userdata');
v = get(btn,'value');
% get data
if ishandle(d.in.handles.hfig)
    D = get(d.in.handles.hfig,'userdata');
    if ~isfield(d.in,'trN')
        trN = 1;
    else
        trN = d.in.trN;
    end
    try
        v = round(v);
        Z = D(d.in.ind,v,trN);
        Z = Z(d.goodChannels);

        clear ud;
        % interpolate data
        ZI = griddata(d.interp.pos(1,:),d.interp.pos(2,:),full(double(Z)),d.interp.XI,d.interp.YI);
        % update data display
        set(d.hi,'Cdata',flipud(ZI));
        % update time index display
        v = round(v);
        set(d.hti,'string',[num2str(d.in.gridTime(v)), ' (', d.in.unit, ')']);
        % update display marker position
        try;set(d.in.hl,'xdata',[v;v]);end
        set(d.ParentAxes,'nextPlot','add')
        try
            % delete current contour plot
            delete(findobj(d.ParentAxes,'type','hggroup'));
            delete(findobj(d.ParentAxes,'type','contour')); % R2014b
            % create new one
            [C,hc] = contour(d.ParentAxes,flipud(ZI),...
                'linecolor',[0.5.*ones(3,1)]);
        end
        axis(d.ParentAxes,'image')
        drawnow
    catch
%     else
        error('Did not find the data!')
    end
else
    error('SPM Graphics Figure has been deleted!')
end
end

%==========================================================================
% get2Dfrom3D
%==========================================================================
function [xy] = get2Dfrom3D(xyz)
% function [xy] = get2Dfrom3D(xyz)
% This function is used to flatten 3D sensor positions onto the 2D plane
% using a modified spherical projection operation.
% It is used to visualize channel data.
% IN:
%   - xyz: the cartesian sensor position in 3D space
% OUT:
%   - xy: the (x,y) cartesian coordinates of the sensors after projection
%   onto the best-fitting sphere

if size(xyz,2) ~= 3
    xyz = xyz';
end
% exclude channels ?
badChannels  = find(isnan(xyz(:,1)));
goodChannels = find(isnan(xyz(:,1))~=1);
xyz          = xyz(goodChannels,:);
% Fit sphere to 3d sensors and center frame
[C,~,~]            = fitSphere(xyz(:,1),xyz(:,2),xyz(:,3));
xyz          = xyz - repmat(C,size(xyz,1),1);
% apply transformation using spherical coordinates
[TH,PHI,RAD] = cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));
TH           = TH - mean(TH);
[X,Y,Z]      = sph2cart(TH,zeros(size(TH)),RAD.*(cos(PHI+pi./2)+1));
xy           = [X(:),Y(:)];
end

%==========================================================================
% fitSphere
%==========================================================================
function [C,R,out] = fitSphere(x,y,z)
% fitSphere  Fit sphere.
%       A = fitSphere(x,y,z) returns the parameters of the best-fit
%       [C,R,out] = fitSphere(x,y,z) returns the center and radius
%       sphere to data points in vectors (x,y,z) using Taubin's method.
% IN:
%   - x/y/z: 3D cartesian coordinates
% OUT:
%   - C: the center of sphere coordinates
%   - R: the radius of the sphere
%   - out: an output structure devoted to graphical display of the best fit
%   sphere

% Make sugary one and zero vectors
l = ones(length(x),1);
O = zeros(length(x),1);

% Make design mx
D = [(x.*x + y.*y + z.*z) x y z l];

Dx = [2*x l O O O];
Dy = [2*y O l O O];
Dz = [2*z O O l O];

% Create scatter matrices
M = D'*D;
N = Dx'*Dx + Dy'*Dy + Dz'*Dz;

% Extract eigensystem
[v, evalues] = eig(M);
evalues = diag(evalues);
Mrank = sum(evalues > eps*5*norm(M));

if (Mrank == 5)
    % Full rank -- min ev corresponds to solution
    Minverse = v'*diag(1./evalues)*v;
    [v,evalues] = eig(inv(M)*N);
    [dmin,dminindex] = max(diag(evalues));
    pvec = v(:,dminindex(1))';
else
    % Rank deficient -- just extract nullspace of M
    pvec = null(M)';
    [m,n] = size(pvec);
    if m > 1
        pvec = pvec(1,:);
    end
end

% Convert to (R,C)
if nargout == 1
    if pvec(1) < 0
        pvec = -pvec;
    end
    C = pvec;
else
    C = -0.5*pvec(2:4) / pvec(1);
    R = sqrt(sum(C*C') - pvec(5)/pvec(1));
end

[X,Y,Z] = sphere;
[TH,PHI,R] = cart2sph(X,Y,Z);
[X,Y,Z] = sph2cart(TH,PHI,R);
X = X + C(1);
Y = Y + C(2);
Z = Z + C(3);

out.X = X;
out.Y = Y;
out.Z = Z;
end
