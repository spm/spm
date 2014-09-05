function [ZI,f] = spm_eeg_plotScalpData(Z,pos,ChanLabel,in)
% Display M/EEG interpolated sensor data on a scalp image
% FORMAT [ZI,f] = spm_eeg_plotScalpData(Z,pos,ChanLabel,in)
%
% INPUT:
%   Z          - the data matrix at the sensors
%   pos        - the positions of the sensors
%   ChanLabel  - the names of the sensors
%   in         - a structure containing some informations related to the 
%                main PRESELECTDATA window. This entry is not necessary
% OUTPUT:
%   ZI         - an image of interpolated data onto the scalp
%   f          - the handle of the figure which displays the interpolated
%                data
%__________________________________________________________________________
%
% This function creates a figure whose purpose is to display an
% interpolation of the sensor data on the scalp (as an image).
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_plotScalpData.m 6157 2014-09-05 18:17:54Z guillaume $


ChanLabel  = char(ChanLabel);
ParentAxes = [];
f          = [];
clim       = [min(Z(:))-( max(Z(:))-min(Z(:)) )/63 , max(Z(:))];
figName    = 'Image Scalp data';
noButtons  = 0;
if nargin < 4 || isempty(in)
    in     = [];
else
    if isfield(in,'min') && ...
            isfield(in,'max') && ...
            isfield(in,'type')
        clim    = [in.min, in.max];
        dc      = abs(diff(clim))./63;
        clim(1) = clim(1) - dc;
        figName = ['Image Scalp data: ',in.type,' sensors'];
        if isfield(in,'trN')
            figName = [figName ', trial #',num2str(in.trN),'.'];
        end
    end
    if isfield(in,'f')
        f = in.f;
    else
        f = figure;
    end
    if isfield(in,'ParentAxes')
        ParentAxes = in.ParentAxes;
    else
        ParentAxes = axes('parent',f);
    end
    if isfield(in,'noButtons')
        noButtons = ~~in.noButtons;
    end     
end

if ~isfield(in,'cbar')
    in.cbar = 1;
end

if ~isfield(in,'plotpos')
    in.plotpos = 1;
end

if size(pos,2) ~= size(ChanLabel, 1)
    pos = pos';
end

nD = size(pos,1);
if nD ~= 2
    % get 2D positions from 3D positions
   xyz   = pos;
   [pos] = get2Dfrom3D(xyz);
   pos   = pos';
end

% exclude channels ?
goodChannels = find(~isnan(pos(1,:)));
pos          = pos(:,goodChannels);
Z            = Z(goodChannels,:);
ChanLabel    = ChanLabel(goodChannels, :);

if ~isempty(in) && isfield(in,'type') && strcmp(in.type, 'MEGPLANAR')
    [cZ, cpos, cChanLabel] = combineplanar(Z, pos, ChanLabel);
else
    cZ         = Z;
    cpos       = pos;
    cChanLabel = ChanLabel;
end

xmin    = min(cpos(1,:));
xmax    = max(cpos(1,:));
dx      = (xmax-xmin)./100;
ymin    = min(cpos(2,:));
ymax    = max(cpos(2,:));
dy      = (ymax-ymin)./100;
x       = xmin:dx:xmax;
y       = ymin:dy:ymax;
[XI,YI] = meshgrid(x,y);
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
fpos(2,:) = 100-fpos(2,:);  % for display purposes (flipud imagesc)

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
    ud = get(in.handles.hfig,'userdata');
    nT = ud.nsamples;
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
if ~noButtons
    set(d.hsp,'userdata',d);
    set(d.hsn,'userdata',d);
end
set(d.ParentAxes,'userdata',d);


%==========================================================================
% dFcn
%==========================================================================
function dFcn(btn,evd)
hf = findobj('tag','Graphics');
D = get(hf,'userdata');
try delete(D.PSD.handles.hli); end


%==========================================================================
% dosp
%==========================================================================
function dosp(btn,evd)
d = get(btn,'userdata');
switch get(d.hp,'visible');
    case 'on'
        set(d.hp,'visible','off');
    case 'off'
        set(d.hp,'visible','on');
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
        Z = D(d.in.ind,v,trN);
        Z = Z(d.goodChannels);

        if strcmp(d.in.type, 'MEGPLANAR')
            Z = combineplanar(Z, d.origpos, d.origChanLabel);
        end

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
C            = fitSphere(xyz(:,1),xyz(:,2),xyz(:,3));
xyz          = xyz - repmat(C,size(xyz,1),1);
% apply transformation using spherical coordinates
[TH,PHI,RAD] = cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));
TH           = TH - mean(TH);
[X,Y,Z]      = sph2cart(TH,zeros(size(TH)),RAD.*(cos(PHI+pi./2)+1));
xy           = [X(:),Y(:)];


%==========================================================================
% combineplanar
%==========================================================================
function [Z, pos, ChanLabel] = combineplanar(Z, pos, ChanLabel)

if ~iscell(ChanLabel)
    ChanLabel = cellstr(ChanLabel);
end

chanind = zeros(1, numel(ChanLabel));
for i = 1:numel(ChanLabel)
    chanind(i) = sscanf(ChanLabel{i}, 'MEG%d');
end

pairs = [];
unpaired = [];
paired = zeros(length(chanind));
for i = 1:length(chanind)
    if ~paired(i)

        cpair = find(abs(chanind - chanind(i))<2);

        if length(cpair) == 1
            unpaired = [unpaired cpair];
        else
            pairs = [pairs; cpair(:)'];
        end
        paired(cpair) = 1;
    end
end

if ~isempty(unpaired)
    warning(['Could not pair all channels. Ignoring ' num2str(length(unpaired)) ' unpaired channels.']);
end

Z = sqrt(Z(pairs(:, 1)).^2 + Z(pairs(:, 2)).^2);
pos = (pos(:, pairs(:, 1)) + pos(:, pairs(:, 2)))./2;
ChanLabel = {};
for i = 1:size(pairs,1)
    ChanLabel{i} = ['MEG' num2str(min(pairs(i,:))) '+' num2str(max(pairs(i,:)))];
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
%   - x/y/z: 3D carthesian ccordinates
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
[TH,PHI,R0] = cart2sph(X,Y,Z);
[X,Y,Z] = sph2cart(TH,PHI,R);
X = X + C(1);
Y = Y + C(2);
Z = Z + C(3);

out.X = X;
out.Y = Y;
out.Z = Z;
